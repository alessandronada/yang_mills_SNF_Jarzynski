#ifndef YM_LOCAL_PT_C
#define YM_LOCAL_PT_C

#include"../include/macro.h"

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>

#ifdef OPENMP_MODE
  #include<omp.h>
#endif

#include"../include/function_pointers.h"
#include"../include/gauge_conf.h"
#include"../include/geometry.h"
#include"../include/gparam.h"
#include"../include/random.h"
#include"../include/spline.h"

double beta_protocol(double t, double b0, double bt)
{
    return b0 + t * (bt - b0);
}

void training_setup(int nsteps, double** forw_plaq, double** forw_work, double** forw_beta)
{
  int err = 0;

  if (err == 0) err = posix_memalign((void**) forw_plaq, DOUBLE_ALIGN, (size_t) (nsteps + 1) * sizeof(double));
  if (err == 0) err = posix_memalign((void**) forw_work, DOUBLE_ALIGN, (size_t) nsteps * sizeof(double));
  if (err == 0) err = posix_memalign((void**) forw_beta, DOUBLE_ALIGN, (size_t) (nsteps + 1) * sizeof(double));

  if (err != 0) {
    fprintf(stderr, "Unable to setup for training %d, %d %s", err, __LINE__, __FILE__);
  }
}

void training_init(int nsteps, double* forw_plaq, double* forw_work)
{
  for (int step = 0; step < nsteps; step++){
    forw_plaq[step] = 0;
    forw_work[step] = 0;
  }
  forw_plaq[nsteps] = 0;
}

void training_clean(double* forw_plaq, double* forw_work, double* forw_beta)
{
  free(forw_work);
  free(forw_plaq);
  free(forw_beta);
}

void real_main(char *in_file)
{
    Gauge_Conf GC, GCstart;
    Geometry geo;
    GParam param;
    double W = 0.0, beta0 = 0.0, dbeta = 0.0, act = 0.0, beta_old = 0.0, plaqs, plaqt;

    // char name[STD_STRING_LENGTH], aux[STD_STRING_LENGTH];
    int count, rel, step, epoch;
    FILE *workfilep;
    time_t time1, time2;

    // to disable nested parallelism
    #ifdef OPENMP_MODE
      // omp_set_nested(0); // deprecated
	omp_set_max_active_levels(1); // should do the same as the old omp_set_nested(0)
    #endif

    // read input file
    readinput(in_file, &param);

    // initialize random generator
    initrand(param.d_randseed);

    // open data_file
    init_work_file(&workfilep, &param);

    // initialize geometry
    init_indexing_lexeo();
    init_geometry(&geo, &param);

    // initialize gauge configuration
    init_gauge_conf(&GC, &param);
    init_gauge_conf_from_gauge_conf(&GCstart, &GC, &param);

    // Monte Carlo begin
    time(&time1);
    beta0 = param.d_beta;
    dbeta = (param.d_J_beta_target - param.d_beta) / param.d_J_steps;

    // thermalization
    for (count = 0; count < param.d_thermal; count++)
    {
        update(&GC, &geo, &param);
    }

    double *forw_plaq, *forw_work, *forw_beta;
    // allocating vectors for training
    training_setup(param.d_J_steps, &forw_plaq, &forw_work, &forw_beta);
    spline* plaq_spline = new_spline(param.d_J_steps+1, forw_beta, forw_plaq);

    forw_beta[0] = beta0;

    // loop on training steps
    for (epoch = 0; epoch < param.d_J_nepochs; epoch++)
    {
        // inint vectors for training
        training_init(param.d_J_steps, forw_plaq, forw_work);
     
        for (step = 0; step < param.d_J_steps; step++) {
          forw_beta[step+1] = beta_protocol(((double) step+1) / ((double) param.d_J_steps), beta0, param.d_J_beta_target);
        }

        // loop on evolutions
        for (count = 0; count < param.d_J_evolutions; count++)
        {
            W = 0.0;
            param.d_beta = beta0;

            // updates between the start of each evolution
            for (rel = 0; rel < param.d_J_between; rel++)
                update(&GC, &geo, &param);

            // increase the index of evolutions
            GC.evolution_index++;

            // store the starting configuration of the evolution
            copy_gauge_conf_from_gauge_conf(&GCstart, &GC, &param);

            // non-equilibrium evolution
            for (step = 0; step < param.d_J_steps; step++)
            {
                //change beta and compute work
                beta_old = param.d_beta;
                param.d_beta = beta_protocol(((double) step+1) / ((double) param.d_J_steps), beta0, param.d_J_beta_target);
                dbeta = param.d_beta - beta_old;
                plaquette(&GC, &geo, &param, &plaqs, &plaqt);
                act = 1 - 0.5 * (plaqs + plaqt);
                act *= (double) (6 * param.d_volume);
                W += dbeta * act;

                // save 0.5 *( plaqs + plaqt ) in array
                forw_plaq[step] += 0.5 * (plaqs + plaqt);

                // saving work increment, full work is computed later
                forw_work[step] += dbeta * act;

                // perform a single step of updates with new beta
                update(&GC, &geo, &param);
            }

            plaquette(&GC, &geo, &param, &plaqs, &plaqt);
            forw_plaq[param.d_J_steps] = 0.5 * (plaqs + plaqt);

            // recover the starting configuration of the evolution
            copy_gauge_conf_from_gauge_conf(&GC, &GCstart, &param);
        }

        double inv_nevols = 1. / param.d_J_evolutions;
        // average of plaquette expectation values for each step
        for (step = 0; step < param.d_J_steps; step++){
          forw_plaq[step] *= inv_nevols;
          forw_work[step] *= inv_nevols;
        }
        forw_plaq[param.d_J_steps] *= inv_nevols;

        // accumulate work increments over steps
        for(step = 1; step < param.d_J_steps; step++) forw_work[step] += forw_work[step - 1];

        // splines interpolation of plaquette expectation values along steps -> P as a function of beta
        init_spline(plaq_spline, param.d_J_steps+1, NULL, NULL);

        // print measures on workfile
        fprintf(workfilep, "%d\n", epoch);
        for(step = 0; step < param.d_J_steps; step++){
          fprintf(workfilep, "%.12g ", forw_beta[step]);
          fprintf(workfilep, "%.12g %.12g ", forw_plaq[step], forw_work[step]);
          fprintf(workfilep, "%.12g\n", integral_spline(plaq_spline, forw_beta[step+1]));
        }
        fprintf(workfilep, "%.12g %.12g\n", forw_beta[param.d_J_steps], forw_plaq[param.d_J_steps]);

        // computation of the gradient of W 

        // protocol parameters update

    }
    // cleaning vectors for training
    free_spline(plaq_spline);
    training_clean(forw_work, forw_plaq, forw_beta);

    time(&time2);
    // Monte Carlo end

    // close data file
    fclose(workfilep);

    // save last starting configuration
    if (param.d_saveconf_back_every != 0)
    {
        write_conf_on_file(&GC, &param);
    }

    // print simulation details
    param.d_beta = beta0;
    print_parameters_j_train_beta(&param, time1, time2);

    // free gauge configurations
    free_gauge_conf(&GC, &param);
    free_gauge_conf(&GCstart, &param);

    // free geometry
    free_geometry(&geo, &param);
}


void print_template_input(void)
  {
  FILE *fp;

  fp=fopen("template_input.example", "w");

  if(fp==NULL)
    {
    fprintf(stderr, "Error in opening the file template_input.example (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  else
    {
    fprintf(fp,"size 4 4 4 4  # Nt Nx Ny Nz\n");
    fprintf(fp,"\n");
    fprintf(fp, "# out-of-equilibrium evolutions parameters\n");
    fprintf(fp, "num_jar_epochs   1         #number of epochs for the training\n");
    fprintf(fp, "num_jar_ev      10         #number of non-equilibrium evolutions\n");
    fprintf(fp, "num_jar_between  1         #number of updates between the start of each evolution\n");
    fprintf(fp, "num_jar_steps   10         #steps in each out-of-equilibrium evolution\n");
    fprintf(fp, "jar_beta_target     6.2    #target beta (only for evolutions in beta)\n");
    fprintf(fp,"\n");
    fprintf(fp,"# Simulations parameters\n");
    fprintf(fp, "beta  5.705\n");
    fprintf(fp, "theta 1.5\n");
    fprintf(fp,"\n");
    fprintf(fp, "thermal    0\n");
    fprintf(fp, "overrelax  5\n");
    fprintf(fp,"\n");
    fprintf(fp, "start                    0  # 0=ordered  1=random  2=from saved configuration\n");
    fprintf(fp, "saveconf_back_every      5  # if 0 does not save, else save backup configurations every ... updates\n");
    fprintf(fp, "saveconf_analysis_every  5  # if 0 does not save, else save configurations for analysis every ... updates\n");
    // fprintf(fp, "\n");
    // fprintf(fp, "coolsteps             3  # number of cooling steps to be used\n");
    // fprintf(fp, "coolrepeat            5  # number of times 'coolsteps' are repeated\n");
    // fprintf(fp, "chi_prime_meas        0  # 1=YES, 0=NO\n");
    // fprintf(fp, "topcharge_tprof_meas  0  # 1=YES, 0=NO\n");
    fprintf(fp,"\n");
    fprintf(fp, "# output files\n");
    fprintf(fp, "conf_file             conf.dat\n");
    fprintf(fp, "log_file              log.dat\n");
    fprintf(fp, "work_file             work.dat\n");
    fprintf(fp, "\n");
    fprintf(fp, "randseed 0    # (0=time)\n");
    fclose(fp);
    }
  }


int main (int argc, char **argv)
    {
    char in_file[STD_STRING_LENGTH];

    if(argc != 2)
      {
			printf("\nSU(N) Hasenbusch Parallel Tempering implemented by Claudio Bonanno (claudiobonanno93@gmail.com) within yang-mills package\n");
			printf("Usage: %s input_file\n\n", argv[0]);

			printf("\nDetails about yang-mills package:\n");
      printf("\tPackage %s version: %s\n", PACKAGE_NAME, PACKAGE_VERSION);
      printf("\tAuthor: Claudio Bonati %s\n\n", PACKAGE_BUGREPORT);

      printf("Compilation details:\n");
      printf("\tN_c (number of colors): %d\n", NCOLOR);
      printf("\tST_dim (space-time dimensionality): %d\n", STDIM);
      printf("\tNum_levels (number of levels): %d\n", NLEVELS);
      printf("\n");
      printf("\tINT_ALIGN: %s\n", QUOTEME(INT_ALIGN));
      printf("\tDOUBLE_ALIGN: %s\n", QUOTEME(DOUBLE_ALIGN));

      #ifdef DEBUG
        printf("\n\tDEBUG mode\n");
      #endif

      #ifdef OPENMP_MODE
        printf("\n\tusing OpenMP with %d threads\n", NTHREADS);
      #endif

      #ifdef THETA_MODE
        printf("\n\tusing imaginary theta\n");
      #endif

      printf("\n");

      #ifdef __INTEL_COMPILER
        printf("\tcompiled with icc\n");
      #elif defined(__clang__)
        printf("\tcompiled with clang\n");
      #elif defined( __GNUC__ )
        printf("\tcompiled with gcc version: %d.%d.%d\n",
                __GNUC__, __GNUC_MINOR__, __GNUC_PATCHLEVEL__);
      #endif

      print_template_input();

      return EXIT_SUCCESS;
      }
    else
      {
      if(strlen(argv[1]) >= STD_STRING_LENGTH)
        {
        fprintf(stderr, "File name too long. Increse STD_STRING_LENGTH in /include/macro.h\n");
				return EXIT_SUCCESS;
        }
      else
        {
    		#if(STDIM==4 && NCOLOR>1)
				strcpy(in_file, argv[1]);
    		real_main(in_file);
    		return EXIT_SUCCESS;
    		#else
    		fprintf(stderr, "Parallel tempering of volume defect not implemented for STDIM =/= 4 and N_color < 2.\n");
    		return EXIT_SUCCESS;
    		#endif
        }
      }
    }

#endif
