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

void real_main(char *in_file)
{
    Gauge_Conf GC, GCstart;
    Geometry geo;
    GParam param;
    double W = 0.0, beta0 = 0.0, old_beta = 0.0, act = 0.0, plaqs, plaqt;

    char name[STD_STRING_LENGTH], aux[STD_STRING_LENGTH];
    int count, rel, step;
    FILE *datafilep, *chiprimefilep,  *topchar_tprof_filep, *workfilep;
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

    // initialize protocol parameters
    init_protocol(&param);

    // open data_file
    init_data_file(&datafilep, &chiprimefilep, &topchar_tprof_filep, &param);
    init_work_file(&workfilep, &param);

    // initialize geometry
    init_indexing_lexeo();
    init_geometry(&geo, &param);

    // initialize gauge configuration
    init_gauge_conf(&GC, &param);
    // copy to save initial configuration on prior
    init_gauge_conf_from_gauge_conf(&GCstart, &GC, &param);

    // Monte Carlo begin
    time(&time1);
    beta0 = param.d_beta;
    //dbeta = (param.d_J_beta_target - param.d_beta) / param.d_J_steps;

    // thermalization
    for (count = 0; count < param.d_thermal; count++)
    {
        update(&GC, &geo, &param);
    }

    // loop on evolutions
    for (count=0; count < param.d_J_evolutions; count++)
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
        old_beta = param.d_beta
        param.d_beta = param.d_J_protocol[step];
          
        plaquette(&GC, &geo, &param, &plaqs, &plaqt);
        act = 1 - 0.5 * (plaqs + plaqt);
        act *= 6 / param.d_inv_vol;
	      W += (param.d_beta - old_beta) * act;
	    
	      // perform a single step of updates with new beta
        update(&GC, &geo, &param);

        if ((step + 1) % param.d_J_dmeas == 0 && step != (param.d_J_steps - 1))
        {
          perform_measures_localobs(&GC, &geo, &param, datafilep, chiprimefilep, topchar_tprof_filep);
          print_work(count, W, workfilep);
        }
	    }

	    // perform measures only on PBC configuration
	    perform_measures_localobs(&GC, &geo, &param, datafilep, chiprimefilep, topchar_tprof_filep);
	    print_work(count, W, workfilep);

      // save initial (beta0) and final (target beta) configurations for offline analysis
      if (param.d_saveconf_analysis_every != 0)
      {
        if (count % param.d_saveconf_analysis_every == 0)
        {
          write_evolution_conf_on_file(&GCstart, &param, 0);
          write_evolution_conf_on_file(&GC, &param, 1);
        }
      }

	    // recover the starting configuration of the evolution
	    copy_gauge_conf_from_gauge_conf(&GC, &GCstart, &param);

      // save initial beta0 configuration for backup
      if (param.d_saveconf_back_every != 0)
      {
        if (count % param.d_saveconf_back_every == 0)
        {
          // simple
          write_conf_on_file(&GC, &param);
          // backup copy
          write_conf_on_file_back(&GC, &param);
        }
      }
    }

    time(&time2);
    // Monte Carlo end

    // close data file
    fclose(datafilep);
    fclose(workfilep);
    if (param.d_chi_prime_meas==1) fclose(chiprimefilep);
    if (param.d_topcharge_tprof_meas==1) fclose(topchar_tprof_filep);

    // save last beta0 configuration
    if (param.d_saveconf_back_every != 0)
    {
      write_conf_on_file(&GC, &param);
    }

    // print simulation details
    param.d_beta = beta0;
    print_parameters_local_jarzynski(&param, time1, time2);

    // free gauge configurations
    free_gauge_conf(&GC, &param);
    free_bound_cond(&GC, &param);
    free_gauge_conf(&GCstart, &param);
    free_bound_cond(&GCstart, &param);

    // free geometry
    free_geometry(&geo, &param);
		
    // free update parameters
    free_hierarc_params(&param);
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
	fprintf(fp,"# parallel tempering parameters\n");
	fprintf(fp,"defect_dir    1             # choose direction of defect boundary: 0->t, 1->x, 2->y, 3->z\n");
	fprintf(fp,"defect_size   1 1 1         # size of the defect (order: y-size z-size t-size)\n");
	fprintf(fp,"N_replica_pt  2    0.0 1.0  # number of parallel tempering replica ____ boundary conditions coefficients\n");
	fprintf(fp,"\n");
    fprintf(fp, "# out-of-equilibrium evolutions parameters\n");
    fprintf(fp, "num_jar_ev      10         #number of non-equilibrium evolutions\n");
    fprintf(fp, "num_jar_between   1        #number of updates between the start of each evolution\n");
    fprintf(fp, "num_jar_steps   10         #steps in each out-of-equilibrium evolution\n");
    fprintf(fp, "num_jar_dmeas   10         #steps between measurements during an evolution (only in beta)\n");
    fprintf(fp, "jar_beta_target     6.2    #target beta (only for evolutions in beta)\n");
	fprintf(fp,"# hierarchical update parameters\n");
	fprintf(fp,"# Ord:qer: num of hierarc levels ____ extension of rectangles ____ num of sweeps per rectangle\n");
	fprintf(fp,"hierarc_upd 2    2 1    1 1\n");
	fprintf(fp,"\n");
	fprintf(fp,"# Simulations parameters\n");
    fprintf(fp, "beta  5.705\n");
    fprintf(fp, "theta 1.5\n");
    fprintf(fp,"\n");
    fprintf(fp, "sample     10\n");
    fprintf(fp, "thermal    0\n");
    fprintf(fp, "overrelax  5\n");
    fprintf(fp, "measevery  1\n");
    fprintf(fp,"\n");
    fprintf(fp, "start                    0  # 0=ordered  1=random  2=from saved configuration\n");
    fprintf(fp, "saveconf_back_every      5  # if 0 does not save, else save backup configurations every ... updates\n");
    fprintf(fp, "saveconf_analysis_every  5  # if 0 does not save, else save configurations for analysis every ... updates\n");
    fprintf(fp, "\n");
	fprintf(fp, "coolsteps             3  # number of cooling steps to be used\n");
	fprintf(fp, "coolrepeat            5  # number of times 'coolsteps' are repeated\n");
	fprintf(fp, "chi_prime_meas        0  # 1=YES, 0=NO\n");
	fprintf(fp, "topcharge_tprof_meas  0  # 1=YES, 0=NO\n");
	fprintf(fp,"\n");
    fprintf(fp, "# output files\n");
	fprintf(fp, "conf_file             conf.dat\n");
	fprintf(fp, "data_file             dati.dat\n");
	fprintf(fp, "chiprime_data_file    chiprime_cool.dat\n");
	fprintf(fp, "topcharge_tprof_file  topo_tcorr_cool.dat\n");
	fprintf(fp, "log_file              log.dat\n");
	fprintf(fp, "swap_acc_file         swap_acc.dat\n");
	fprintf(fp, "swap_track_file       swap_track.dat\n");
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
