#ifndef YM_LOCAL_SNF_BETA_C
#define YM_LOCAL_SNF_BETA_C

#include "../include/macro.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#ifdef OPENMP_MODE
#include <omp.h>
#endif

#include "../include/function_pointers.h"
#include "../include/gauge_conf.h"
#include "../include/geometry.h"
#include "../include/gparam.h"
#include "../include/random.h"

void real_main(char *in_file)
{
  Gauge_Conf GC, GCstart;
  Geometry geo;
  GParam param;
  double W = 0.0, beta_0 = 0.0, beta_t_0 = 0.0, act0 = 0.0, act1 = 0.0, plaqs, plaqt, logJ;

  int npar, count, rel, step;
  FILE *datafilep, *chiprimefilep, *topchar_tprof_filep, *workfilep;
  time_t time1, time2;

  // to disable nested parallelism
  #ifdef OPENMP_MODE
    // omp_set_nested(0); // deprecated
    omp_set_max_active_levels(1); // should do the same as the old omp_set_nested(0)
  #endif

  // read input file
  readinput(in_file, &param);
  check_flow_beta_input(&param);

  // initialize random generator
  initrand(param.d_randseed);

  // initialize protocol parameters
  if (param.d_anisotropic)
    npar = 2;
  else
    npar = 1;
  init_start_end_protocol_beta(&param, npar);
  init_protocol(&param, npar);

  // initialize smearing parameters
  init_smearing_parameter(&param);

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
  beta_0 = param.d_beta;
  if (param.d_anisotropic != 0)
    beta_t_0 = param.d_beta_t;

  // thermalization
  for (count = 0; count < param.d_thermal; count++)
  {
    update(&GC, &geo, &param);
  }

  // loop on evolutions
  for (count = 0; count < param.d_flow_evolutions; count++)
  {
    W = 0.0;
    param.d_beta = beta_0;
    if (param.d_anisotropic != 0)
      param.d_beta_t = beta_t_0;

    // updates between the start of each evolution
    for (rel = 0; rel < param.d_flow_between; rel++)
      update(&GC, &geo, &param);

    // increase the index of evolutions
    GC.evolution_index++;

    // store the starting configuration of the evolution
    copy_gauge_conf_from_gauge_conf(&GCstart, &GC, &param);

    // non-equilibrium evolution
    for (step = 0; step < param.d_flow_steps; step++)
    {
      // compute S_beta(i) (U_i)
      plaquette(&GC, &geo, &param, &plaqs, &plaqt);
      act0 = wilson_action(&param, param.d_beta, param.d_beta_t, plaqs, plaqt);

      // stout smearing step: U_i -> g_i(U_i)
      isotropic_stout_smearing_update(&GC, &geo, &param, &logJ, (param.d_SNF_rho)[step]);

      // change beta: S_beta(i) -> S_beta(i+1)
      set_flow_beta(&param, step);

      // compute S_beta(i+1) (g_i(U_i)) and work
      plaquette(&GC, &geo, &param, &plaqs, &plaqt);
      act1 = wilson_action(&param, param.d_beta, param.d_beta_t, plaqs, plaqt);
      W += act1 - act0 - logJ;

      // MC update: g_i(U_i) -> U_(i+1)
      update(&GC, &geo, &param);

      if ((step + 1) % param.d_flow_dmeas == 0 && step != (param.d_flow_steps - 1))
      {
        perform_measures_localobs(&GC, &geo, &param, datafilep, chiprimefilep, topchar_tprof_filep);
        print_work((int)GC.evolution_index, W, workfilep);
      }
    }

    // perform measures only on PBC configuration
    perform_measures_localobs(&GC, &geo, &param, datafilep, chiprimefilep, topchar_tprof_filep);
    print_work((int)GC.evolution_index, W, workfilep);

    // save initial (beta0) and final (target beta) configurations for offline analysis
    if (param.d_saveconf_analysis_every != 0)
    {
      if ((int)GC.evolution_index % param.d_saveconf_analysis_every == 0)
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
  if (param.d_chi_prime_meas == 1)
    fclose(chiprimefilep);
  if (param.d_topcharge_tprof_meas == 1)
    fclose(topchar_tprof_filep);

  // save last beta0 configuration
  if (param.d_saveconf_back_every != 0)
  {
    write_conf_on_file(&GC, &param);
  }

  // print simulation details
  param.d_beta = beta_0;
  if (param.d_anisotropic != 0)
    param.d_beta_t = beta_t_0;
  print_parameters_local_flow_beta(&param, time1, time2);

  // free gauge configurations
  free_gauge_conf(&GC, &param);
  free_gauge_conf(&GCstart, &param);

  // free geometry
  free_geometry(&geo, &param);

  // free protocol and smearing parameters
  free_flow_params(&param);
}

void print_template_input(void)
{
  FILE *fp;

  fp = fopen("template_input.example", "w");

  if (fp == NULL)
  {
    fprintf(stderr, "Error in opening the file template_input.example (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  else
  {
    fprintf(fp, "size 4 4 4 4  # Nt Nx Ny Nz\n");
    fprintf(fp, "\n");
    fprintf(fp, "# action (beta, beta_t: couplings at the start of each evolution)\n");
    fprintf(fp, "beta         5.705\n");
    fprintf(fp, "anisotropic  0      # 0 = isotropic, otherwise beta for spatial and beta_t for temporal plaquettes\n");
    fprintf(fp, "beta_t       5.705  # (only if anisotropic)\n");
    fprintf(fp, "theta        0.0    # imaginary theta (only if compiled with --enable-use-theta)\n");
    fprintf(fp, "\n");
    fprintf(fp, "# flow in beta (SNF: stout smearing and update at each step)\n");
    fprintf(fp, "flow_beta_target    6.2  # target beta\n");
    fprintf(fp, "flow_beta_t_target  6.2  # target beta_t (only if anisotropic)\n");
    fprintf(fp, "num_flow_ev         10   # number of non-equilibrium evolutions\n");
    fprintf(fp, "num_flow_between    1    # number of updates between the start of two evolutions\n");
    fprintf(fp, "num_flow_steps      10   # number of steps of each evolution\n");
    fprintf(fp, "num_flow_dmeas      10   # steps between measurements during an evolution\n");
    fprintf(fp, "protocol_type       0    # 0 = linear protocol, otherwise read from protocol_file\n");
    fprintf(fp, "\n");
    fprintf(fp, "# Monte Carlo\n");
    fprintf(fp, "thermal    0\n");
    fprintf(fp, "overrelax  5\n");
    fprintf(fp, "\n");
    fprintf(fp, "start                    0  # 0=ordered  1=random  2=from saved configuration\n");
    fprintf(fp, "saveconf_back_every      5  # if 0 does not save, else save backup configurations every ... evolutions\n");
    fprintf(fp, "saveconf_analysis_every  5  # if 0 does not save, else save configurations for analysis every ... evolutions\n");
    fprintf(fp, "\n");
    fprintf(fp, "# measurements\n");
    fprintf(fp, "coolsteps             3  # number of cooling steps to be used\n");
    fprintf(fp, "coolrepeat            5  # number of times 'coolsteps' are repeated\n");
    fprintf(fp, "chi_prime_meas        0  # 1=YES, 0=NO\n");
    fprintf(fp, "topcharge_tprof_meas  0  # 1=YES, 0=NO\n");
    fprintf(fp, "\n");
    fprintf(fp, "# input and output files\n");
    fprintf(fp, "smearingrho_file      rho.dat              # rho of the stout smearing, one per step (text)\n");
    fprintf(fp, "conf_file             conf.dat\n");
    fprintf(fp, "data_file             dati.dat\n");
    fprintf(fp, "work_file             work.dat\n");
    fprintf(fp, "protocol_file         protocol.dat         # (only if protocol_type != 0)\n");
    fprintf(fp, "chiprime_data_file    chiprime_cool.dat    # (only if chi_prime_meas = 1)\n");
    fprintf(fp, "topcharge_tprof_file  topo_tcorr_cool.dat  # (only if topcharge_tprof_meas = 1)\n");
    fprintf(fp, "log_file              log.dat\n");
    fprintf(fp, "\n");
    fprintf(fp, "randseed 0    # (0=time)\n");
    fclose(fp);
  }
}

int main(int argc, char **argv)
{
  char in_file[STD_STRING_LENGTH];

  if (argc != 2)
  {
    printf("\nNE-MCMC and SNF in beta implemented by Alessandro Nada (nada.alessandro@gmail.com)\n");
    printf("\nStout smearing routines implemented along with Dario Panfalone and Lorenzo Verzichelli\n");
    printf("Usage: %s input_file\n\n", argv[0]);

    printf("\nForked from yang_mills_PTBC by Claudio Bonanno (claudiobonanno93@gmail.com) within yang-mills package\n");

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
#elif defined(__GNUC__)
    printf("\tcompiled with gcc version: %d.%d.%d\n",
           __GNUC__, __GNUC_MINOR__, __GNUC_PATCHLEVEL__);
#endif

    print_template_input();

    return EXIT_SUCCESS;
  }
  else
  {
    if (strlen(argv[1]) >= STD_STRING_LENGTH)
    {
      fprintf(stderr, "File name too long. Increase STD_STRING_LENGTH in /include/macro.h\n");
      return EXIT_FAILURE;
    }
    else
    {
#if (STDIM == 4 && NCOLOR == 3)
      strcpy(in_file, argv[1]);
      real_main(in_file);
      return EXIT_SUCCESS;
#else
      fprintf(stderr, "SNF implemented only for STDIM = 4 and N_c = 3 (Jacobian of the stout smearing).\n");
      return EXIT_FAILURE;
#endif
    }
  }
}

#endif
