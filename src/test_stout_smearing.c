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
#include"../include/sun.h"

void print_SUN(SuN const * const A) {
    for (int i = 0; i < NCOLOR; i++) {
        for (int j = 0; j < NCOLOR; j++) {
            printf("%lf + %lfi ", creal(A->comp[m(i,j)]), cimag(A->comp[m(i,j)]));
        }
        printf("\n");
    }
    printf("\n");
}

void test_taexp_SU3(SuN const * const test_A) {
    printf("Testing taexp for SU(3)\n");
    SuN SUN_expA, SU3_expA;

    // print_SUN(test_A);

    equal_SuN(&SUN_expA, test_A);
    taexp_SuN(&SUN_expA);
    // print_SUN(&SUN_expA);

    equal_SuN(&SU3_expA, test_A);
    taexp(&SU3_expA); // should call taexp_Su3
    // print_SUN(&SU3_expA);

    double elem_max = 0.;
    for (int i = 0; i < NCOLOR * NCOLOR; i++) {
        // printf("%lf + %lfi, %lf + %lfi \n",
        //     creal(SUN_expA.comp[i]), cimag(SUN_expA.comp[i]), creal(SU3_expA.comp[i]), cimag(SU3_expA.comp[i]));
        double elem = cabs(SUN_expA.comp[i] - SU3_expA.comp[i]);
        if (elem != elem) printf("NaN enecountered! \n");
        if (elem > elem_max) {
            elem_max = elem;
        }
    }

    printf("Max difference between SUN_expA and SU3_expA: %e\n", elem_max);
}

void isotropic_smear_conf_for_test(Gauge_Conf const * const source,
                              Gauge_Conf* smeared,
                              GParam const * const param,
                              Geometry const * const geo,
                              double const rho) {
    smeared->update_index = source->update_index;
    for (int dir = 0; dir < STDIM; dir++) {
        for (long r = 0; r < (param->d_volume)/2; r++) {
            isotropic_stout_smearing_singlelink(source, geo, param, r, dir, rho, &(smeared->lattice[r][dir]));
        }
        for (long r = (param->d_volume)/2; r < param->d_volume; r++) {
            isotropic_stout_smearing_singlelink(source, geo, param, r, dir, rho, &(smeared->lattice[r][dir]));
        }
    }
}

void anisotropic_smear_conf_for_test(Gauge_Conf const * const source,
                                     Gauge_Conf* smeared,
                                     GParam const * const param,
                                     Geometry const * const geo,
                                     double *rho) {
    smeared->update_index = source->update_index;
    for (int dir = 0; dir < STDIM; dir++) {
        for (long r = 0; r < (param->d_volume)/2; r++) {
            anisotropic_stout_smearing_singlelink(source, geo, r, dir, rho + dir * STDIM, &(smeared->lattice[r][dir]));
        }
        for (long r = (param->d_volume)/2; r < param->d_volume; r++) {
            anisotropic_stout_smearing_singlelink(source, geo, r, dir, rho + dir * STDIM, &(smeared->lattice[r][dir]));
        }
    }
}

void real_main(char* in_file) {
    
    GParam param;
    readinput(in_file, &param);
    initrand(param.d_randseed);

    print_parameters_local(&param, 0, 0);

    Geometry geo;
    init_indexing_lexeo();
    init_geometry(&geo, &param);

    Gauge_Conf GC, smeared_GC;
    init_gauge_conf(&GC, &param);
    init_gauge_conf_from_gauge_conf(&smeared_GC, &GC, &param);

    // Test the taexp function for SU(3)
    SuN A;
    one(&A);
    A.comp[1] += 0.1;
    test_taexp_SU3(&A);

    rand_matrix_SuN(&A);
    test_taexp_SU3(&A);

    // Add any additional tests or functionality here
    double rho[STDIM*STDIM];
    for (int mu = 0; mu < STDIM; mu++) {
        printf("%d ", mu); fflush(stdout);
        for (int nu = 0; nu < STDIM; nu++) {
            rho[mu*STDIM + nu] = casuale();
        }
    }
    // printf("%f \n", rho[0]);

    printf("%f ", creal(GC.lattice[0][0].comp[0]));
    isotropic_smear_conf_for_test(&GC, &smeared_GC, &param, &geo, rho[0]);
    printf("%f ", creal(smeared_GC.lattice[0][0].comp[0]));
    anisotropic_smear_conf_for_test(&GC, &smeared_GC, &param, &geo, rho);
    printf("%f \n", creal(smeared_GC.lattice[0][0].comp[0]));

    for (int i = 1; i <= param.d_thermal; i++){
        update(&GC, &geo, &param);

        if (param.d_saveconf_back_every != 0)
            if (i % param.d_saveconf_back_every == 0) {
                write_conf_on_file_back(&GC, &param);
        }
    }

    FILE *datafilep, *chipfilep, *topcfilep;
    init_data_file(&datafilep, &chipfilep, &topcfilep, &param);
    double rho_value = 0.;

    for (int i = 1; i <= param.d_sample; i++) {
        update(&GC, &geo, &param);

        if (param.d_measevery != 0)
            if (i % param.d_measevery == 0) {
                double plaqs, plaqt;
                plaquette(&GC, &geo, &param, &plaqs, &plaqt);
                fprintf(datafilep, "%ld %.3f %e %e \n", GC.update_index, 0., plaqs, plaqt);
                for (rho_value = 0.; rho_value < 0.5; rho_value += 0.01) {
                    isotropic_smear_conf_for_test(&GC, &smeared_GC, &param, &geo, rho_value);
                    plaquette(&smeared_GC, &geo, &param, &plaqs, &plaqt);
                    fprintf(datafilep, "%ld %.3f %e %e \n", GC.update_index, rho_value, plaqs, plaqt);
                }
                fflush(datafilep);
        }

        if (param.d_saveconf_back_every != 0)
            if (i % param.d_saveconf_back_every == 0) {
                write_conf_on_file(&GC, &param);
        }
    }

    fclose(datafilep);
    if (param.d_chi_prime_meas) fclose(chipfilep);
    if (param.d_topcharge_tprof_meas) fclose(topcfilep);

    free_gauge_conf(&GC, &param);
    free_gauge_conf(&smeared_GC, &param);
    free_geometry(&geo, &param);
}

int main(int argc, char **argv) {
    if (NCOLOR != 3) {
        fprintf(stderr, "This test is for SU(3)\n");
        return EXIT_FAILURE;
    }
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

    //   print_template_input();

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
            strcpy(in_file, argv[1]);
    		real_main(in_file);
    		return EXIT_SUCCESS;
        }
      }
    }
