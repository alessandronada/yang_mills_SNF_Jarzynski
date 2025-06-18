#ifndef GPARAM_C
#define GPARAM_C

#include"../include/macro.h"
#include"../include/endianness.h"
#include"../include/gparam.h"

#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>


// remove from input file white/empty lines and comments
// comments start with the charachter #
void remove_white_line_and_comments(FILE *input)
  {
  int temp_i;

  temp_i=getc(input);
  if(temp_i=='\n' || temp_i==' ' || temp_i=='\043') // scan for white lines and comments
    {
    ungetc(temp_i, input);

    temp_i=getc(input);
    if(temp_i=='\n' || temp_i==' ') // white line
      {
      do
       {
       temp_i=getc(input);
       }
      while(temp_i=='\n' || temp_i==' ');
      }
    ungetc(temp_i, input);

    temp_i=getc(input);
    if(temp_i=='\043')  // comment, \043 = ascii oct for #
      {
      do
       {
       temp_i=getc(input);
       }
      while(temp_i!='\n');
      }
    else
      {
      ungetc(temp_i, input);
      }

    remove_white_line_and_comments(input);
    }
  else
    {
    ungetc(temp_i, input);
    }
  }


void readinput(char *in_file, GParam *param)
    {
    FILE *input;
    char str[STD_STRING_LENGTH], temp_str[STD_STRING_LENGTH];
    double temp_d;
    int temp_i, i;
    int err, end=1;
    unsigned int temp_ui;

    // this is to avoid unnecessary checks in case the multilevel is not used
    for(i=0; i<NLEVELS; i++)
       {
       param->d_ml_step[i]=0;
       }

    // just to avoid possible mistakes with uninitialized stuff
    for(i=0; i<NCOLOR; i++)
       {
       param->d_h[i]=0.0;
       }
    param->d_theta=0.0;

		// to avoid possible mistakes with uninitialized stuff
		for (i=0; i<STDIM-1; i++)
		{
		param->d_L_defect[i]=0;
		}
		param->d_N_replica_pt=1;

		param->d_J_evolutions = 0;
		param->d_J_between = 0;
		param->d_J_steps = 0;
    param->d_J_dmeas = 0;
    param->d_J_beta_target = 6.0;
    param->d_J_protocol_type = 0;
		
		// default = do not compute chi_prime
		param->d_chi_prime_meas = 0;
		param->d_topcharge_tprof_meas = 0;

    input=fopen(in_file, "r"); // open the input file
    if(input==NULL)
      {
      fprintf(stderr, "Error in opening the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
      exit(EXIT_FAILURE);
      }
    else
      {
      while(end==1) // slide the file
           {
           remove_white_line_and_comments(input);

           err=fscanf(input, "%s", str);
           if(err!=1)
             {
             fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
             printf("err=%d\n", err);
             exit(EXIT_FAILURE);
             }

           if(strncmp(str, "size", 4)==0)
             {
             for(i=0; i<STDIM; i++)
                {
                err=fscanf(input, "%d", &temp_i);
                if(err!=1)
                  {
                  fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
                  exit(EXIT_FAILURE);
                  }
                param->d_size[i]=temp_i;
                }
             }

           else if(strncmp(str, "beta", 4)==0)
                  { 
                  err=fscanf(input, "%lf", &temp_d);
                  if(err!=1)
                    {
                    fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                    }
                  param->d_beta=temp_d;
                  }
           else if(strncmp(str, "htracedef", 9)==0)
                  {
                  for(i=0; i<(int)floor(NCOLOR/2.0); i++)
                     {
                     err=fscanf(input, "%lf", &temp_d);
                     if(err!=1)
                       {
                       fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
                       exit(EXIT_FAILURE);
                       }
                     param->d_h[i]=temp_d;
                     }
                  }
           else if(strncmp(str, "theta", 5)==0)
                  {
                  err=fscanf(input, "%lf", &temp_d);
                  if(err!=1)
                    {
                    fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                    }
                  param->d_theta=temp_d;
                  }

           else if(strncmp(str, "sample", 6)==0)
                  { 
                  err=fscanf(input, "%d", &temp_i);
                  if(err!=1)
                    {
                    fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                    }
                  param->d_sample=temp_i;
                  }
           else if(strncmp(str, "thermal", 7)==0)
                  { 
                  err=fscanf(input, "%d", &temp_i);
                  if(err!=1)
                    {
                    fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                    }
                  param->d_thermal=temp_i;
                  }
           else if(strncmp(str, "overrelax", 9)==0)
                  { 
                  err=fscanf(input, "%d", &temp_i);
                  if(err!=1)
                    {
                    fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                    }
                  param->d_overrelax=temp_i;
                  }
           else if(strncmp(str, "measevery", 9)==0)
                  { 
                  err=fscanf(input, "%d", &temp_i);
                  if(err!=1)
                    {
                    fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                    }
                  param->d_measevery=temp_i;
                  }

           else if(strncmp(str, "start", 5)==0)
                  { 
                  err=fscanf(input, "%d", &temp_i);
                  if(err!=1)
                    {
                    fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                    }
                  param->d_start=temp_i;
                  }
           else if(strncmp(str, "saveconf_back_every", 19)==0)
                  {
                  err=fscanf(input, "%d", &temp_i);
                  if(err!=1)
                    {
                    fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                    }
                  param->d_saveconf_back_every=temp_i;
                  }
           else if(strncmp(str, "saveconf_analysis_every", 23)==0)
                  {
                  err=fscanf(input, "%d", &temp_i);
                  if(err!=1)
                    {
                    fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                    }
                  param->d_saveconf_analysis_every=temp_i;
                  }

           else if(strncmp(str, "epsilon_metro", 13)==0)
                  {
                  err=fscanf(input, "%lf", &temp_d);
                  if(err!=1)
                    {
                    fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                    }
                  param->d_epsilon_metro=temp_d;
                  }

           else if(strncmp(str, "coolsteps", 9)==0)
                  {
                  err=fscanf(input, "%d", &temp_i);
                  if(err!=1)
                    {
                    fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                    }
                  param->d_coolsteps=temp_i;
                  }
           else if(strncmp(str, "coolrepeat", 10)==0)
                  {
                  err=fscanf(input, "%d", &temp_i);
                  if(err!=1)
                    {
                    fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                    }
                  param->d_coolrepeat=temp_i;
                  }
           else if(strncmp(str, "chi_prime_meas", 14)==0)
                  {
                  err=fscanf(input, "%d", &temp_i);
                  if(err!=1)
                    {
                    fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                    }
									if ( (temp_i == 1) || (temp_i == 0 ) ) param->d_chi_prime_meas=temp_i;
									else
										{
											fprintf(stderr, "Error: chi_prime_meas must be either 0 or 1 in %s (%s, %d)\n", in_file, __FILE__, __LINE__);
											exit(EXIT_FAILURE);
										}
                  }

           else if(strncmp(str, "topcharge_tprof_meas", 20)==0)
                  {
                  err=fscanf(input, "%d", &temp_i);
                  if(err!=1)
                    {
                    fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                    }
									if ( (temp_i == 1) || (temp_i == 0 ) ) param->d_topcharge_tprof_meas=temp_i;
									else
										{
											fprintf(stderr, "Error: chi_prime_meas must be either 0 or 1 in %s (%s, %d)\n", in_file, __FILE__, __LINE__);
											exit(EXIT_FAILURE);
										}
                  }
									
           else if(strncmp(str, "topcharge_tprof_file", 20)==0)
                  { 
                  err=fscanf(input, "%s", temp_str);
                  if(err!=1)
                    {
                    fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                    }
                  strcpy(param->d_topcharge_tprof_file, temp_str);
                  }

           else if(strncmp(str, "gfstep", 6)==0) // integration step
                  {
                  err=fscanf(input, "%lf", &temp_d);
                  if(err!=1)
                    {
                    fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                    }
                  param->d_gfstep=temp_d;
                  }

           else if(strncmp(str, "num_gfsteps", 11)==0) // number of integration steps
                  {
                  err=fscanf(input, "%d", &temp_i);
                  if(err!=1)
                    {
                    fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                    }
                  param->d_ngfsteps=temp_i;
                  }
           else if(strncmp(str, "gf_meas_each", 12)==0) // number of integration steps
                  {
                  err=fscanf(input, "%d", &temp_i);
                  if(err!=1)
                    {
                    fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                    }
                  param->d_gf_meas_each=temp_i;
                  }

           else if(strncmp(str, "multihit", 8)==0)
                  {
                  err=fscanf(input, "%d", &temp_i);
                  if(err!=1)
                    {
                    fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                    }
                  param->d_multihit=temp_i;
                  }
           else if(strncmp(str, "ml_step", 7)==0)
                  {
                  for(i=0; i<NLEVELS; i++)
                     {
                     err=fscanf(input, "%d", &temp_i);
                     if(err!=1)
                       {
                       fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
                       exit(EXIT_FAILURE);
                       }
                     param->d_ml_step[i]=temp_i;
                     }
                  }
           else if(strncmp(str, "ml_upd", 6)==0)
                  {
                  for(i=0; i<NLEVELS; i++)
                     {
                     err=fscanf(input, "%d", &temp_i);
                     if(err!=1)
                       {
                       fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
                       exit(EXIT_FAILURE);
                       }
                     param->d_ml_upd[i]=temp_i;
                     }
                  }
           else if(strncmp(str, "ml_level0_repeat", 16)==0)
                  {
                  err=fscanf(input, "%d", &temp_i);
                  if(err!=1)
                    {
                    fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                    }
                  param->d_ml_level0_repeat=temp_i;
                  }
           else if(strncmp(str, "dist_poly", 9)==0)
                  {
                  err=fscanf(input, "%d", &temp_i);
                  if(err!=1)
                    {
                    fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                    }
                  param->d_dist_poly=temp_i;
                  }
           else if(strncmp(str, "transv_dist", 11)==0)
                  {
                  err=fscanf(input, "%d", &temp_i);
                  if(err!=1)
                    {
                    fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                    }
                  param->d_trasv_dist=temp_i;
                  }
           else if(strncmp(str, "plaq_dir", 8)==0)
                  {
                  err=fscanf(input, "%d", &temp_i);
                  if(err!=1)
                    {
                    fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                    }
                  param->d_plaq_dir[0]=temp_i;

                  err=fscanf(input, "%d", &temp_i);
                  if(err!=1)
                    {
                    fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                    }
                  param->d_plaq_dir[1]=temp_i;
                  }

           else if(strncmp(str, "conf_file", 9)==0)
                  { 
                  err=fscanf(input, "%s", temp_str);
                  if(err!=1)
                    {
                    fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                    }
                  strcpy(param->d_conf_file, temp_str);
                  }
           else if(strncmp(str, "data_file", 9)==0)
                  { 
                  err=fscanf(input, "%s", temp_str);
                  if(err!=1)
                    {
                    fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                    }
                  strcpy(param->d_data_file, temp_str);
                  }
	    else if(strncmp(str, "work_file", 9)==0)
                  {
                  err=fscanf(input, "%s", temp_str);
                  if(err!=1)
                    {
                    fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                    }
                  strcpy(param->d_work_file, temp_str);
                  }
      else if(strncmp(str, "protocol_file", 13)==0)
                  {
                  err=fscanf(input, "%s", temp_str);
                  if(err!=1)
                    {
                    fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                    }
                  strcpy(param->d_protocol_file, temp_str);
                  }
      else if(strncmp(str, "smearingrho_file", 16)==0)
                  {
                  err=fscanf(input, "%s", temp_str);
                  if(err!=1)
                    {
                    fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                    }
                  strcpy(param->d_smearingrho_file, temp_str);
                  }
           else if(strncmp(str, "chiprime_data_file", 18)==0)
                  { 
                  err=fscanf(input, "%s", temp_str);
                  if(err!=1)
                    {
                    fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                    }
                  strcpy(param->d_chiprime_file, temp_str);
                  }
           else if(strncmp(str, "log_file", 8)==0)
                  { 
                  err=fscanf(input, "%s", temp_str);
                  if(err!=1)
                    {
                    fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                    }
                  strcpy(param->d_log_file, temp_str);
                  }
           else if(strncmp(str, "ml_file", 7)==0)
                  {
                  err=fscanf(input, "%s", temp_str);
                  if(err!=1)
                    {
                    fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                    }
                  strcpy(param->d_ml_file, temp_str);
                  }

           else if(strncmp(str, "randseed", 8)==0)
                  { 
                  err=fscanf(input, "%u", &temp_ui);
                  if(err!=1)
                    {
                    fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                    }
                  param->d_randseed=temp_ui;
                  }
				else if(strncmp(str, "defect_dir", 10)==0)
                  {
                  err=fscanf(input, "%d", &temp_i);
                  if(err!=1)
                    {
                    fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                    }
                  if( (temp_i==0) || (temp_i==1) || (temp_i==2) || (temp_i==3) )
					{
					param->d_defect_dir=temp_i;
					}
				  else
					{
					fprintf(stderr, "Error in reading the file %s, defect_dir must be either 0,1,2 or 3 (%s, %d)\n", in_file, __FILE__, __LINE__);
					exit(EXIT_FAILURE);
					}
                  }
						else if(strncmp(str, "defect_size", 11)==0)
									{
									for (i=0; i<STDIM-1; i++)
									{
									err=fscanf(input, "%d", &temp_i);
									if(err!=1)
										{
										fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
										exit(EXIT_FAILURE);
										}
									param->d_L_defect[i]=temp_i;
									}
									}
           else if(strncmp(str, "N_replica_pt", 12)==0)
                  {
                  err=fscanf(input, "%d", &temp_i);
                  if(err!=1)
                    {
                    fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                    }
									param->d_N_replica_pt=temp_i;
									err=posix_memalign( (void **) &(param->d_pt_bound_cond_coeff), (size_t) DOUBLE_ALIGN, (size_t) param->d_N_replica_pt * sizeof(double));
									if(err!=0)
										{
										fprintf(stderr, "Problems in allocating parallel tempering parameters! (%s, %d)\n", __FILE__, __LINE__);
										exit(EXIT_FAILURE);
										}
									for(i=0;i<param->d_N_replica_pt;i++)
										{
										err=fscanf(input, "%lf", &temp_d);
										if(err!=1)
											{
											fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
											exit(EXIT_FAILURE);
											}
										param->d_pt_bound_cond_coeff[i]=temp_d;
										}
                  }
           else if (strncmp(str, "jar_beta_target", 15) == 0)
           {
                  err = fscanf(input, "%lf", &temp_d);
                  if (err != 1)
                  {
                      fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
                      exit(EXIT_FAILURE);
                  }
                  param->d_J_beta_target = temp_d;
           }
		   else if (strncmp(str, "num_jar_ev", 10) == 0)
		   {
				err = fscanf(input, "%d", &temp_i);
				if (err != 1)
				{
					fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
					exit(EXIT_FAILURE);
				}
				param->d_J_evolutions = temp_i;
		   }
		   else if (strncmp(str, "num_jar_steps", 13) == 0)
		   {
				err = fscanf(input, "%d", &temp_i);
				if (err != 1)
				{
					fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
					exit(EXIT_FAILURE);
				}
			param->d_J_steps = temp_i;
		   }
		   else if (strncmp(str, "num_jar_between", 13) == 0)
			{
			err = fscanf(input, "%d", &temp_i);
			if (err != 1)
			{
				fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
				exit(EXIT_FAILURE);
			}
			param->d_J_between = temp_i;
			}
           else if (strncmp(str, "num_jar_dmeas", 13) == 0)
            {
            err = fscanf(input, "%d", &temp_i);
            if (err != 1)
            {
                fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
                exit(EXIT_FAILURE);
            }
            param->d_J_dmeas = temp_i;
            }
            else if (strncmp(str, "protocol_type", 13) == 0)
            {
              err = fscanf(input, "%d", &temp_i);
              if (err != 1)
              {
                fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
                exit(EXIT_FAILURE);
              }
              param->d_J_protocol_type = temp_i;
            }
           else if(strncmp(str, "swap_acc_file", 13)==0)
                  { 
                  err=fscanf(input, "%s", temp_str);
                  if(err!=1)
                    {
                    fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                    }
                  strcpy(param->d_swap_acc_file, temp_str);
                  }			  
           else if(strncmp(str, "swap_track_file", 15)==0)
                  { 
                  err=fscanf(input, "%s", temp_str);
                  if(err!=1)
                    {
                    fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                    }
                  strcpy(param->d_swap_tracking_file, temp_str);
                  }
           else if(strncmp(str, "hierarc_upd", 11)==0)
                  { 
                  err=fscanf(input, "%d", &temp_i);
                  if(err!=1)
                    {
                    fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                    }
                  param->d_N_hierarc_levels=temp_i;
									if(param->d_N_hierarc_levels > 0)
										{
										err=posix_memalign( (void **) &(param->d_L_rect), (size_t) INT_ALIGN, (size_t) param->d_N_hierarc_levels * sizeof(int));
										if (err!=0)
											{
											fprintf(stderr, "Problems in allocating hierarchical update parameters! (%s, %d)\n", __FILE__, __LINE__);
											exit(EXIT_FAILURE);
											}
										for(i=0;i<param->d_N_hierarc_levels;i++)
											{
											err=fscanf(input, "%d", &temp_i);
											if(err!=1)
												{
												fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
												exit(EXIT_FAILURE);
												}
											param->d_L_rect[i]=temp_i;
											}
										err=posix_memalign( (void **) &(param->d_N_sweep_rect), (size_t) INT_ALIGN, (size_t) param->d_N_hierarc_levels * sizeof(int));
										if (err!=0)
											{
											fprintf(stderr, "Problems in allocating hierarchical update parameters! (%s, %d)\n", __FILE__, __LINE__);
											exit(EXIT_FAILURE);
											}
										for(i=0;i<param->d_N_hierarc_levels;i++)
											{
											err=fscanf(input, "%d", &temp_i);
											if(err!=1)
												{
												fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
												exit(EXIT_FAILURE);
												}
											param->d_N_sweep_rect[i]=temp_i;
											}
										} // closes if( num_hierarc_levels > 0 )
                  }
							else if(strncmp(str, "multicanonic_acc_file", 21)==0)
							{ 
								err=fscanf(input, "%s", temp_str);
								if(err!=1)
								{
								fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
								exit(EXIT_FAILURE);
								}
								strcpy(param->d_multicanonic_acc_file, temp_str);
							}
							else if(strncmp(str, "topo_potential_file", 19)==0)
							{ 
								err=fscanf(input, "%s", temp_str);
								if(err!=1)
								{
									fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
									exit(EXIT_FAILURE);
								}
								strcpy(param->d_topo_potential_file, temp_str);
							}
							else if(strncmp(str, "grid_step", 9)==0)
							{
								err=fscanf(input, "%lf", &temp_d);
								if(err!=1)
								{
									fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
									exit(EXIT_FAILURE);
								}
								param->d_grid_step=temp_d;
							}
							else if(strncmp(str, "grid_max", 8)==0)
							{
								err=fscanf(input, "%lf", &temp_d);
								if(err!=1)
								{
									fprintf(stderr, "Error in reading the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
									exit(EXIT_FAILURE);
								}
								param->d_grid_max=temp_d;
							}
           else
             {
             fprintf(stderr, "Error: unrecognized option %s in the file %s (%s, %d)\n", str, in_file, __FILE__, __LINE__);
             exit(EXIT_FAILURE);
             }

           remove_white_line_and_comments(input);

           // check if the read line is the last one
           temp_i=getc(input);
           if(temp_i==EOF)
             {
             end=0;
             }
           else
             {
             ungetc(temp_i, input);
             }
           }

      fclose(input);

      // VARIOUS CHECKS
      if(param->d_ml_step[0]!=0)
        {
        if(param->d_size[0] % param->d_ml_step[0] || param->d_size[0] < param->d_ml_step[0])
          {
          fprintf(stderr, "Error: size[0] has to be divisible by ml_step[0] and satisfy ml_step[0]<=size[0] (%s, %d)\n", __FILE__, __LINE__);
          exit(EXIT_FAILURE);
          }
        for(i=1; i<NLEVELS; i++)
           {
           if(param->d_ml_step[i-1] % param->d_ml_step[i] || param->d_ml_step[i-1] <= param->d_ml_step[i])
             {
             fprintf(stderr, "Error: ml_step[%d] has to be divisible by ml_step[%d] and larger than it (%s, %d)\n", i-1, i, __FILE__, __LINE__);
             exit(EXIT_FAILURE);
             }
           }
        if(param->d_ml_step[NLEVELS-1]==1)
          {
          fprintf(stderr, "Error: ml_step[%d] has to be larger than 1 (%s, %d)\n", NLEVELS-1, __FILE__, __LINE__);
          exit(EXIT_FAILURE);
          }
        }

      #ifdef OPENMP_MODE
      for(i=0; i<STDIM; i++)
         {
         temp_i = param->d_size[i] % 2;
         if(temp_i!=0)
           {
           fprintf(stderr, "Error: size[%d] is not even.\n", i);
           fprintf(stderr, "When using OpenMP all the sides of the lattice have to be even! (%s, %d)\n", __FILE__, __LINE__);
           exit(EXIT_FAILURE);
           }
         }
      #endif

      err=0;
      for(i=0; i<STDIM; i++)
         {
         if(param->d_size[i]==1)
           {
           err=1;
           }
         }
      if(err==1)
        {
        fprintf(stderr, "Error: all sizes has to be larger than 1: the totally reduced case is not implemented! (%s, %d)\n", __FILE__, __LINE__);
        }
				
			// various checks on parallel tempering parameters
			if(param->d_L_defect[0]>param->d_size[0])
				{
				fprintf(stderr, "Error: defect's t-length is greater than lattice's t-length (%s, %d)\n", __FILE__, __LINE__);
				exit(EXIT_FAILURE);
				}
			if(param->d_L_defect[1]>param->d_size[2])
				{
				fprintf(stderr, "Error: defect's y-length is greater than lattice's y-length (%s, %d)\n", __FILE__, __LINE__);
				exit(EXIT_FAILURE);
				}
			if(param->d_L_defect[2]>param->d_size[3])
				{
				fprintf(stderr, "Error: defect's z-length is greater than lattice's z-length (%s, %d)\n", __FILE__, __LINE__);
				exit(EXIT_FAILURE);
				}
			if(param->d_N_replica_pt<1)
				{
				fprintf(stderr, "Error: number of replica of parallel tempering must be greater than 0 (%s, %d)\n", __FILE__, __LINE__);
				exit(EXIT_FAILURE);
				}

      init_derived_constants(param);
      }
    }


void init_derived_constants(GParam *param)
  {
  int i;

  // derived constants
  param->d_volume=1;
  for(i=0; i<STDIM; i++)
     {
     (param->d_volume)*=(param->d_size[i]);
     }

  param->d_space_vol=1;
  // direction 0 is time
  for(i=1; i<STDIM; i++)
     {
     (param->d_space_vol)*=(param->d_size[i]);
     }

  param->d_inv_vol=1.0/((double) param->d_volume);
  param->d_inv_space_vol=1.0/((double) param->d_space_vol);
  
  // volume of the defect
  param->d_volume_defect=1;
  for(i=0; i<STDIM-1;i++)
  {
	param->d_volume_defect *= param->d_L_defect[i];
  }
	
	// number of grid points (multicanonic only)
	param->d_n_grid=(int)((2.0*param->d_grid_max/param->d_grid_step)+1.0);
  }

  

void init_protocol(GParam const * const param, double start, double end)
{
  FILE *input_protocol;
  double temp_d;
  int i;
  int err;

  err=posix_memalign( (void **) &(param->d_J_protocol), (size_t) DOUBLE_ALIGN, (size_t) param->d_J_steps * sizeof(double));
	if(err!=0)
	{
	  fprintf(stderr, "Problems in allocating protocol parameters! (%s, %d)\n", __FILE__, __LINE__);
	  exit(EXIT_FAILURE);
	}			

  if (param->d_J_protocol_type)
  {
    input_protocol=fopen(param->d_protocol_file, "r"); // open the input protocol file
    
    if(input_protocol==NULL)
    {
      fprintf(stderr, "Error in opening the file %s (%s, %d)\n", param->d_protocol_file, __FILE__, __LINE__);
      exit(EXIT_FAILURE);
    }
    else
    {
      for(i=0;i<param->d_J_steps;i++)
	    {
        err=fscanf(input_protocol, "%lf", &temp_d);
			  if(err!=1)
			  {
			    fprintf(stderr, "Error in reading the file %s (%s, %d)\n", param->d_protocol_file, __FILE__, __LINE__);
			    exit(EXIT_FAILURE);
			  }
			  param->d_J_protocol[i]=temp_d;
      }
    }
  }
  else
  {
    for(i=0;i<param->d_J_steps;i++)
      param->d_J_protocol[i]=(double)((end - start) * ((double)(i+1)) / param->d_J_steps + start);
  }
}

void init_smearing_parameter(GParam const * const param)
{
  FILE *input_smearingrho;
  double temp_d;
  int i;
  int err;

  err=posix_memalign( (void **) &(param->d_SNF_rho), (size_t) DOUBLE_ALIGN, (size_t) param->d_J_steps * sizeof(double));
	if(err!=0)
	{
	  fprintf(stderr, "Problems in allocating protocol parameters! (%s, %d)\n", __FILE__, __LINE__);
	  exit(EXIT_FAILURE);
	}			

  input_smearingrho=fopen(param->d_smearingrho_file, "r"); // open the input smearing rho file

  if(input_smearingrho==NULL)
  {
    fprintf(stderr, "Error in opening the file %s (%s, %d)\n", param->d_smearingrho_file, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  else
  {
    for(i=0;i<param->d_J_steps;i++)
	  {
      err=fscanf(input_smearingrho, "%lf", &temp_d);
	    if(err!=1)
		  {
		    fprintf(stderr, "Error in reading the file %s (%s, %d)\n", param->d_smearingrho_file, __FILE__, __LINE__);
		    exit(EXIT_FAILURE);
		  }
		  param->d_SNF_rho[i]=temp_d;

    }
  }
}

void init_defect_smearing_parameter(GParam const * const param, long rect_vol)
{
  FILE *input_smearingrho;
  double temp_d;
  int i, mu, s, p;
  int err;

  err=posix_memalign( (void **) &(param->d_SNF_rho), (size_t) DOUBLE_ALIGN, (size_t) param->d_J_steps * 2 * (STDIM-1) * STDIM * rect_vol * sizeof(double));
	if(err!=0)
	{
	  fprintf(stderr, "Problems in allocating protocol parameters! (%s, %d)\n", __FILE__, __LINE__);
	  exit(EXIT_FAILURE);
	}			

  input_smearingrho=fopen(param->d_smearingrho_file, "r"); // open the input smearing rho file

  if(input_smearingrho==NULL)
  {
    fprintf(stderr, "Error in opening the file %s (%s, %d)\n", param->d_smearingrho_file, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  else
  {
    for(i=0;i<param->d_J_steps;i++)
      for(mu=0;mu<STDIM;mu++)
        for(s=0;s<rect_vol;s++)
          for(p=0;p<2*(STDIM-1);p++)
          {
            long rho_index = 2 * (STDIM-1) * rect_vol * STDIM * i + 2 * (STDIM-1) * rect_vol * mu + 2 * (STDIM-1) * s + p;
            err=fscanf(input_smearingrho, "%lf", &temp_d);
	          if(err!=1)
            { 
              fprintf(stderr, "Error in reading the file %s (%s, %d)\n", param->d_smearingrho_file, __FILE__, __LINE__);
              exit(EXIT_FAILURE);
            }
		        param->d_SNF_rho[rho_index]=temp_d;
          }
  }
}

// initialize data file
void init_data_file(FILE **dataf, FILE **chiprimef, FILE **topchar_tprof_f, GParam const * const param)
{
	int i;

	if(param->d_start==2)
	{
		// open std data file (plaquette, polyakov, topological charge)
		*dataf=fopen(param->d_data_file, "r");
		if(*dataf!=NULL) // file exists
		{
			fclose(*dataf);
			*dataf=fopen(param->d_data_file, "a");
		}
		else
		{
 			*dataf=fopen(param->d_data_file, "w");
			fprintf(*dataf, "%d ", STDIM);
			for(i=0; i<STDIM; i++) fprintf(*dataf, "%d ", param->d_size[i]);
      fprintf(*dataf, "\n");
		}
		// open chi prime data file
		if (param->d_chi_prime_meas == 1)
		{
			*chiprimef=fopen(param->d_chiprime_file, "r");
			if(*chiprimef!=NULL) // file exists
			{
				fclose(*chiprimef);
				*chiprimef=fopen(param->d_chiprime_file, "a");
			}
			else
			{
 				*chiprimef=fopen(param->d_chiprime_file, "w");
				fprintf(*chiprimef, "%d ", STDIM);
				for(i=0; i<STDIM; i++) fprintf(*chiprimef, "%d ", param->d_size[i]);
      	fprintf(*chiprimef, "\n");
			}
		}
		else
		{
			(void) chiprimef;
		}
		// open topocharge_tprof data file
		if (param->d_topcharge_tprof_meas == 1)
		{
			*topchar_tprof_f=fopen(param->d_topcharge_tprof_file, "r");
			if(*topchar_tprof_f!=NULL) // file exists
			{
				fclose(*topchar_tprof_f);
				*topchar_tprof_f=fopen(param->d_topcharge_tprof_file, "a");
			}
			else
			{
 				*topchar_tprof_f=fopen(param->d_topcharge_tprof_file, "w");
				fprintf(*topchar_tprof_f, "%d ", STDIM);
				for(i=0; i<STDIM; i++) fprintf(*topchar_tprof_f, "%d ", param->d_size[i]);
      	fprintf(*topchar_tprof_f, "\n");
			}
		}
		else
		{
			(void) topchar_tprof_f;
		}
	}
	else
	{
		// open std data file
    *dataf=fopen(param->d_data_file, "w");
    fprintf(*dataf, "%d ", STDIM);
    for(i=0; i<STDIM; i++)
       {
       fprintf(*dataf, "%d ", param->d_size[i]);
       }
    fprintf(*dataf, "\n");
		// open chi prime data file
		if (param->d_chi_prime_meas == 1)
		{
    	*chiprimef=fopen(param->d_chiprime_file, "w");
    	fprintf(*chiprimef, "%d ", STDIM);
    	for(i=0; i<STDIM; i++) fprintf(*chiprimef, "%d ", param->d_size[i]);
    	fprintf(*chiprimef, "\n");
		}
		else
		{
			(void) chiprimef;
		}
		// open topocharge_tprof data file
		if (param->d_topcharge_tprof_meas == 1)
		{
    	*topchar_tprof_f=fopen(param->d_topcharge_tprof_file, "w");
    	fprintf(*topchar_tprof_f, "%d ", STDIM);
    	for(i=0; i<STDIM; i++) fprintf(*topchar_tprof_f, "%d ", param->d_size[i]);
    	fprintf(*topchar_tprof_f, "\n");
		}
		else
		{
			(void) topchar_tprof_f;
		}
	}
	fflush(*dataf);
	if (param->d_chi_prime_meas == 1 ) fflush(*chiprimef);
	else 
		{
			(void) chiprimef;
		}
	if (param->d_topcharge_tprof_meas == 1 ) fflush(*topchar_tprof_f);
	else 
		{
			(void) topchar_tprof_f;
		}
}

void init_work_file(FILE **workfilep, GParam const * const param)
{
	*workfilep = fopen(param->d_work_file, "r");
	if (*workfilep != NULL) // file exists
	{
		fclose(*workfilep);
		*workfilep = fopen(param->d_work_file, "a");
	}
	else
	{
		int i;
		*workfilep = fopen(param->d_work_file, "w");
		fprintf(*workfilep, "# %f ", param->d_beta);
		fprintf(*workfilep, "%d ", STDIM);
		for (i = 0; i < STDIM; i++) fprintf(*workfilep, "%d ", param->d_size[i]);
		fprintf(*workfilep, "\n");
	}
	fflush(*workfilep);
}

// free allocated memory for hierarc update parameters
void free_hierarc_params(GParam *param)
	{
	if(param->d_N_hierarc_levels==0)
		{
		(void) param; // to avoid compiler warning about unused variable
		}
	else
		{
		free(param->d_L_rect);
		free(param->d_N_sweep_rect);
		}
	}


// print simulation parameters

void print_parameters_local(GParam const * const param, time_t time_start, time_t time_end)
    {
    FILE *fp;
    int i;
    double diff_sec;

    fp=fopen(param->d_log_file, "w");
    fprintf(fp, "+-----------------------------------------+\n");
    fprintf(fp, "| Simulation details for yang_mills_local |\n");
    fprintf(fp, "+-----------------------------------------+\n\n");

    #ifdef OPENMP_MODE
     fprintf(fp, "using OpenMP with %d threads\n\n", NTHREADS);
    #endif

    fprintf(fp, "number of colors: %d\n", NCOLOR);
    fprintf(fp, "spacetime dimensionality: %d\n\n", STDIM);

    fprintf(fp, "lattice: %d", param->d_size[0]);
    for(i=1; i<STDIM; i++)
       {
       fprintf(fp, "x%d", param->d_size[i]);
       }
    fprintf(fp, "\n\n");

    fprintf(fp, "beta: %.10lf\n", param->d_beta);
    #ifdef THETA_MODE
      fprintf(fp, "theta: %.10lf\n", param->d_theta);
    #endif
    fprintf(fp, "\n");

    fprintf(fp, "sample:    %d\n", param->d_sample);
    fprintf(fp, "thermal:   %d\n", param->d_thermal);
    fprintf(fp, "overrelax: %d\n", param->d_overrelax);
    fprintf(fp, "measevery: %d\n", param->d_measevery);
    fprintf(fp, "\n");

    fprintf(fp, "start:                   %d\n", param->d_start);
    fprintf(fp, "saveconf_back_every:     %d\n", param->d_saveconf_back_every);
    fprintf(fp, "saveconf_analysis_every: %d\n", param->d_saveconf_analysis_every);
    fprintf(fp, "\n");

    fprintf(fp, "coolsteps:      %d\n", param->d_coolsteps);
    fprintf(fp, "coolrepeat:     %d\n", param->d_coolrepeat);
    fprintf(fp, "\n");

    fprintf(fp, "randseed: %u\n", param->d_randseed);
    fprintf(fp, "\n");

    diff_sec = difftime(time_end, time_start);
    fprintf(fp, "Simulation time: %.3lf seconds\n", diff_sec );
    fprintf(fp, "\n");

    if(endian()==0)
      {
      fprintf(fp, "Little endian machine\n\n");
      }
    else
      {
      fprintf(fp, "Big endian machine\n\n");
      }

    fclose(fp);
    }

void print_parameters_local_pt_multicanonic(GParam const * const param, time_t time_start, time_t time_end)
    {
    FILE *fp;
    int i;
    double diff_sec;

    fp=fopen(param->d_log_file, "w");
    fprintf(fp, "+---------------------------------------------------------+\n");
    fprintf(fp, "| Simulation details for yang_mills_local_pt_multicanonic |\n");
    fprintf(fp, "+---------------------------------------------------------+\n\n");

    #ifdef OPENMP_MODE
     fprintf(fp, "using OpenMP with %d threads\n\n", NTHREADS);
    #endif

    fprintf(fp, "number of colors: %d\n", NCOLOR);
    fprintf(fp, "spacetime dimensionality: %d\n\n", STDIM);

    fprintf(fp, "lattice: %d", param->d_size[0]);
    for(i=1; i<STDIM; i++)
       {
       fprintf(fp, "x%d", param->d_size[i]);
       }
    fprintf(fp, "\n\n");
	
	fprintf(fp, "defect dir: %d\n", param->d_defect_dir);
    fprintf(fp, "defect: %d", param->d_L_defect[0]);
    for(i=1; i<STDIM-1; i++)
       {
       fprintf(fp, "x%d", param->d_L_defect[i]);
       }
    fprintf(fp, "\n\n");
    fprintf(fp,"number of copies used in parallel tempering: %d\n", param->d_N_replica_pt);
		fprintf(fp,"boundary condition constants: ");
		for(i=0;i<param->d_N_replica_pt;i++)
			fprintf(fp,"%lf ",param->d_pt_bound_cond_coeff[i]);
		fprintf(fp,"\n");
		fprintf(fp,"number of hierarchical levels: %d\n", param->d_N_hierarc_levels);
		if(param->d_N_hierarc_levels>0)
			{
			fprintf(fp,"extention of rectangles: ");
			for(i=0;i<param->d_N_hierarc_levels;i++)
				{
				fprintf(fp,"%d ", param->d_L_rect[i]);
				}
			fprintf(fp,"\n");
			fprintf(fp,"number of sweeps per hierarchical level: ");
			for(i=0;i<param->d_N_hierarc_levels;i++)
				{
				fprintf(fp,"%d ", param->d_N_sweep_rect[i]);
				}
			}
		fprintf(fp,"\n\n");
		
		fprintf(fp,"Multicanonic topo-potential read from file %s\nPotential defined on a grid with step=%.10lf and max=%.10lf\n", param->d_topo_potential_file, param->d_grid_step, param->d_grid_max);
		
		fprintf(fp,"\n\n");

    fprintf(fp, "beta: %.10lf\n", param->d_beta);
    #ifdef THETA_MODE
      fprintf(fp, "theta: %.10lf\n", param->d_theta);
    #endif
    fprintf(fp, "\n");

    fprintf(fp, "sample:    %d\n", param->d_sample);
    fprintf(fp, "thermal:   %d\n", param->d_thermal);
    fprintf(fp, "overrelax: %d\n", param->d_overrelax);
    fprintf(fp, "measevery: %d\n", param->d_measevery);
    fprintf(fp, "\n");

    fprintf(fp, "start:                   %d\n", param->d_start);
    fprintf(fp, "saveconf_back_every:     %d\n", param->d_saveconf_back_every);
    fprintf(fp, "saveconf_analysis_every: %d\n", param->d_saveconf_analysis_every);
    fprintf(fp, "\n");

    fprintf(fp, "coolsteps:      %d\n", param->d_coolsteps);
    fprintf(fp, "coolrepeat:     %d\n", param->d_coolrepeat);
    fprintf(fp, "\n");

    fprintf(fp, "randseed: %u\n", param->d_randseed);
    fprintf(fp, "\n");

    diff_sec = difftime(time_end, time_start);
    fprintf(fp, "Simulation time: %.3lf seconds\n", diff_sec );
    fprintf(fp, "\n");

    if(endian()==0)
      {
      fprintf(fp, "Little endian machine\n\n");
      }
    else
      {
      fprintf(fp, "Big endian machine\n\n");
      }

    fclose(fp);
    }

void print_parameters_local_pt(GParam const * const param, time_t time_start, time_t time_end)
    {
    FILE *fp;
    int i;
    double diff_sec;

    fp=fopen(param->d_log_file, "w");
    fprintf(fp, "+--------------------------------------------+\n");
    fprintf(fp, "| Simulation details for yang_mills_local_pt |\n");
    fprintf(fp, "+--------------------------------------------+\n\n");

    #ifdef OPENMP_MODE
     fprintf(fp, "using OpenMP with %d threads\n\n", NTHREADS);
    #endif

    fprintf(fp, "number of colors: %d\n", NCOLOR);
    fprintf(fp, "spacetime dimensionality: %d\n\n", STDIM);

    fprintf(fp, "lattice: %d", param->d_size[0]);
    for(i=1; i<STDIM; i++)
       {
       fprintf(fp, "x%d", param->d_size[i]);
       }
    fprintf(fp, "\n\n");
	
	fprintf(fp, "defect dir: %d\n", param->d_defect_dir);
    fprintf(fp, "defect: %d", param->d_L_defect[0]);
    for(i=1; i<STDIM-1; i++)
       {
       fprintf(fp, "x%d", param->d_L_defect[i]);
       }
    fprintf(fp, "\n\n");
    fprintf(fp,"number of copies used in parallel tempering: %d\n", param->d_N_replica_pt);
		fprintf(fp,"boundary condition constants: ");
		for(i=0;i<param->d_N_replica_pt;i++)
			fprintf(fp,"%lf ",param->d_pt_bound_cond_coeff[i]);
		fprintf(fp,"\n");
		fprintf(fp,"number of hierarchical levels: %d\n", param->d_N_hierarc_levels);
		if(param->d_N_hierarc_levels>0)
			{
			fprintf(fp,"extention of rectangles: ");
			for(i=0;i<param->d_N_hierarc_levels;i++)
				{
				fprintf(fp,"%d ", param->d_L_rect[i]);
				}
			fprintf(fp,"\n");
			fprintf(fp,"number of sweeps per hierarchical level: ");
			for(i=0;i<param->d_N_hierarc_levels;i++)
				{
				fprintf(fp,"%d ", param->d_N_sweep_rect[i]);
				}
			}
		fprintf(fp,"\n\n");

    fprintf(fp, "beta: %.10lf\n", param->d_beta);
    #ifdef THETA_MODE
      fprintf(fp, "theta: %.10lf\n", param->d_theta);
    #endif
    fprintf(fp, "\n");

    fprintf(fp, "sample:    %d\n", param->d_sample);
    fprintf(fp, "thermal:   %d\n", param->d_thermal);
    fprintf(fp, "overrelax: %d\n", param->d_overrelax);
    fprintf(fp, "measevery: %d\n", param->d_measevery);
    fprintf(fp, "\n");

    fprintf(fp, "start:                   %d\n", param->d_start);
    fprintf(fp, "saveconf_back_every:     %d\n", param->d_saveconf_back_every);
    fprintf(fp, "saveconf_analysis_every: %d\n", param->d_saveconf_analysis_every);
    fprintf(fp, "\n");

    fprintf(fp, "coolsteps:      %d\n", param->d_coolsteps);
    fprintf(fp, "coolrepeat:     %d\n", param->d_coolrepeat);
    fprintf(fp, "\n");

    fprintf(fp, "randseed: %u\n", param->d_randseed);
    fprintf(fp, "\n");

    diff_sec = difftime(time_end, time_start);
    fprintf(fp, "Simulation time: %.3lf seconds\n", diff_sec );
    fprintf(fp, "\n");

    if(endian()==0)
      {
      fprintf(fp, "Little endian machine\n\n");
      }
    else
      {
      fprintf(fp, "Big endian machine\n\n");
      }

    fclose(fp);
    }

void print_parameters_local_jarzynski(GParam const * const param, time_t time_start, time_t time_end)
{
	FILE *fp;
	int i;
	double diff_sec;

	fp = fopen(param->d_log_file, "w");
	fprintf(fp, "+---------------------------------------------------+\n");
	fprintf(fp, "| Simulation details for yang_mills_local_jarzynski |\n");
	fprintf(fp, "+---------------------------------------------------+\n\n");

#ifdef OPENMP_MODE
	fprintf(fp, "using OpenMP with %d threads\n\n", NTHREADS);
#endif

	fprintf(fp, "number of colors: %d\n", NCOLOR);
	fprintf(fp, "spacetime dimensionality: %d\n\n", STDIM);

	fprintf(fp, "lattice: %d", param->d_size[0]);
	for (i = 1; i < STDIM; i++)
	{
		fprintf(fp, "x%d", param->d_size[i]);
	}
	fprintf(fp, "\n\n");

	fprintf(fp, "defect dir: %d\n", param->d_defect_dir);
	fprintf(fp, "defect: %d", param->d_L_defect[0]);
	for (i = 1; i < STDIM - 1; i++)
	{
		fprintf(fp, "x%d", param->d_L_defect[i]);
	}
	fprintf(fp, "\n\n");
	fprintf(fp, "number of out-of-equilibrium evolutions: %d\n", param->d_J_evolutions);
	fprintf(fp, "number of out-of-equilibrium steps in each evolution: %d\n", param->d_J_steps);
	fprintf(fp, "number of relax updates between evolutions: %d\n", param->d_J_between);
    fprintf(fp, "number of steps between measurements during evolution: %d\n", param->d_J_dmeas);
	fprintf(fp, "number of hierarchical levels: %d\n", param->d_N_hierarc_levels);
	if (param->d_N_hierarc_levels > 0)
	{
		fprintf(fp, "extention of rectangles: ");
		for (i = 0; i < param->d_N_hierarc_levels; i++)
		{
			fprintf(fp, "%d ", param->d_L_rect[i]);
		}
		fprintf(fp, "\n");
		fprintf(fp, "number of sweeps per hierarchical level: ");
		for (i = 0; i < param->d_N_hierarc_levels; i++)
		{
			fprintf(fp, "%d ", param->d_N_sweep_rect[i]);
		}
	}
	fprintf(fp, "\n\n");

	fprintf(fp, "beta: %.10lf\n", param->d_beta);
    fprintf(fp, "beta target: %.10lf\n", param->d_J_beta_target);
#ifdef THETA_MODE
	fprintf(fp, "theta: %.10lf\n", param->d_theta);
#endif
	fprintf(fp, "\n");

	fprintf(fp, "sample:    %d\n", param->d_sample);
	fprintf(fp, "thermal:   %d\n", param->d_thermal);
	fprintf(fp, "overrelax: %d\n", param->d_overrelax);
	fprintf(fp, "measevery: %d\n", param->d_measevery);
	fprintf(fp, "\n");

	fprintf(fp, "start:                   %d\n", param->d_start);
	fprintf(fp, "saveconf_back_every:     %d\n", param->d_saveconf_back_every);
	fprintf(fp, "saveconf_analysis_every: %d\n", param->d_saveconf_analysis_every);
	fprintf(fp, "\n");

	fprintf(fp, "coolsteps:      %d\n", param->d_coolsteps);
	fprintf(fp, "coolrepeat:     %d\n", param->d_coolrepeat);
	fprintf(fp, "\n");

	fprintf(fp, "randseed: %u\n", param->d_randseed);
	fprintf(fp, "\n");

	diff_sec = difftime(time_end, time_start);
	fprintf(fp, "Simulation time: %.3lf seconds\n", diff_sec);
	fprintf(fp, "\n");

	if (endian() == 0)
	{
		fprintf(fp, "Little endian machine\n\n");
	}
	else
	{
		fprintf(fp, "Big endian machine\n\n");
	}

	fclose(fp);
}

// print simulation parameters
void print_parameters_polycorr_long(GParam * param, time_t time_start, time_t time_end)
    {
    FILE *fp;
    int i;
    double diff_sec;

    fp=fopen(param->d_log_file, "w");
    fprintf(fp, "+-------------------------------------------------+\n");
    fprintf(fp, "| Simulation details for yang_mills_polycorr_long |\n");
    fprintf(fp, "+-------------------------------------------------+\n\n");

    #ifdef OPENMP_MODE
     fprintf(fp, "using OpenMP with %d threads\n\n", NTHREADS);
    #endif

    fprintf(fp, "number of colors: %d\n", NCOLOR);
    fprintf(fp, "spacetime dimensionality: %d\n\n", STDIM);

    fprintf(fp, "lattice: %d", param->d_size[0]);
    for(i=1; i<STDIM; i++)
       {
       fprintf(fp, "x%d", param->d_size[i]);
       }
    fprintf(fp, "\n\n");

    fprintf(fp, "beta: %.10lf\n", param->d_beta);
    #ifdef THETA_MODE
      fprintf(fp, "theta: %.10lf\n", param->d_theta);
    #endif
    fprintf(fp, "\n");

    fprintf(fp, "sample:    %d\n", param->d_sample);
    fprintf(fp, "thermal:   %d\n", param->d_thermal);
    fprintf(fp, "overrelax: %d\n", param->d_overrelax);
    fprintf(fp, "measevery: %d\n", param->d_measevery);
    fprintf(fp, "\n");

    fprintf(fp, "start:                   %d\n", param->d_start);
    fprintf(fp, "\n");

    fprintf(fp, "multihit:   %d\n", param->d_multihit);
    fprintf(fp, "levels for multileves: %d\n", NLEVELS);
    fprintf(fp, "multilevel steps: ");
    for(i=0; i<NLEVELS; i++)
       {
       fprintf(fp, "%d ", param->d_ml_step[i]);
       }
    fprintf(fp, "\n");
    fprintf(fp, "updates for levels: ");
    for(i=0; i<NLEVELS; i++)
       {
       fprintf(fp, "%d ", param->d_ml_upd[i]);
       }
    fprintf(fp, "\n");
    fprintf(fp, "level0_repeat:   %d\n", param->d_ml_level0_repeat);
    fprintf(fp, "dist_poly:   %d\n", param->d_dist_poly);
    fprintf(fp, "\n");

    fprintf(fp, "randseed: %u\n", param->d_randseed);
    fprintf(fp, "\n");

    diff_sec = difftime(time_end, time_start);
    fprintf(fp, "Simulation time: %.3lf seconds\n", diff_sec );
    fprintf(fp, "\n");

    if(endian()==0)
      {
      fprintf(fp, "Little endian machine\n\n");
      }
    else
      {
      fprintf(fp, "Big endian machine\n\n");
      }

    fclose(fp);
    }


// print simulation parameters
void print_parameters_polycorr(GParam * param, time_t time_start, time_t time_end)
    {
    FILE *fp;
    int i;
    double diff_sec;

    fp=fopen(param->d_log_file, "w");
    fprintf(fp, "+--------------------------------------------+\n");
    fprintf(fp, "| Simulation details for yang_mills_polycorr |\n");
    fprintf(fp, "+--------------------------------------------+\n\n");

    #ifdef OPENMP_MODE
     fprintf(fp, "using OpenMP with %d threads\n\n", NTHREADS);
    #endif

    fprintf(fp, "number of colors: %d\n", NCOLOR);
    fprintf(fp, "spacetime dimensionality: %d\n\n", STDIM);

    fprintf(fp, "lattice: %d", param->d_size[0]);
    for(i=1; i<STDIM; i++)
       {
       fprintf(fp, "x%d", param->d_size[i]);
       }
    fprintf(fp, "\n\n");

    fprintf(fp, "beta: %.10lf\n", param->d_beta);
    #ifdef THETA_MODE
      fprintf(fp, "theta: %.10lf\n", param->d_theta);
    #endif

    fprintf(fp, "\n");

    fprintf(fp, "sample:    %d\n", param->d_sample);
    fprintf(fp, "thermal:   %d\n", param->d_thermal);
    fprintf(fp, "overrelax: %d\n", param->d_overrelax);
    fprintf(fp, "measevery: %d\n", param->d_measevery);
    fprintf(fp, "\n");

    fprintf(fp, "start:                   %d\n", param->d_start);
    fprintf(fp, "saveconf_back_every:     %d\n", param->d_saveconf_back_every);
    fprintf(fp, "\n");

    fprintf(fp, "multihit:   %d\n", param->d_multihit);
    fprintf(fp, "levels for multileves: %d\n", NLEVELS);
    fprintf(fp, "multilevel steps: ");
    for(i=0; i<NLEVELS; i++)
       {
       fprintf(fp, "%d ", param->d_ml_step[i]);
       }
    fprintf(fp, "\n");
    fprintf(fp, "updates for levels: ");
    for(i=0; i<NLEVELS; i++)
       {
       fprintf(fp, "%d ", param->d_ml_upd[i]);
       }
    fprintf(fp, "\n");
    fprintf(fp, "dist_poly:  %d\n", param->d_dist_poly);
    fprintf(fp, "\n");

    fprintf(fp, "randseed: %u\n", param->d_randseed);
    fprintf(fp, "\n");

    diff_sec = difftime(time_end, time_start);
    fprintf(fp, "Simulation time: %.3lf seconds\n", diff_sec );
    fprintf(fp, "\n");

    if(endian()==0)
      {
      fprintf(fp, "Little endian machine\n\n");
      }
    else
      {
      fprintf(fp, "Big endian machine\n\n");
      }

    fclose(fp);
    }


// print simulation parameters
void print_parameters_t0(GParam * param, time_t time_start, time_t time_end)
    {
    FILE *fp;
    int i;
    double diff_sec;

    fp=fopen(param->d_log_file, "w");
    fprintf(fp, "+--------------------------------------+\n");
    fprintf(fp, "| Simulation details for yang_mills_t0 |\n");
    fprintf(fp, "+--------------------------------------+\n\n");

    #ifdef OPENMP_MODE
     fprintf(fp, "using OpenMP with %d threads\n\n", NTHREADS);
    #endif

    fprintf(fp, "number of colors: %d\n", NCOLOR);
    fprintf(fp, "spacetime dimensionality: %d\n\n", STDIM);

    fprintf(fp, "lattice: %d", param->d_size[0]);
    for(i=1; i<STDIM; i++)
       {
       fprintf(fp, "x%d", param->d_size[i]);
       }
    fprintf(fp, "\n\n");

    fprintf(fp, "gfstep:    %lf\n", param->d_gfstep);
    fprintf(fp, "\n");

    fprintf(fp, "randseed: %u\n", param->d_randseed);
    fprintf(fp, "\n");

    diff_sec = difftime(time_end, time_start);
    fprintf(fp, "Simulation time: %.3lf seconds\n", diff_sec );
    fprintf(fp, "\n");

    if(endian()==0)
      {
      fprintf(fp, "Little endian machine\n\n");
      }
    else
      {
      fprintf(fp, "Big endian machine\n\n");
      }

    fclose(fp);
    }

// print simulation parameters
void print_parameters_gf(GParam * param, time_t time_start, time_t time_end)
    {
    FILE *fp;
    int i;
    double diff_sec;

    fp=fopen(param->d_log_file, "w");
    fprintf(fp, "+-------------------------------------------------+\n");
    fprintf(fp, "| Simulation details for yang_mills_gradient_flow |\n");
    fprintf(fp, "+-------------------------------------------------+\n\n");

    #ifdef OPENMP_MODE
     fprintf(fp, "using OpenMP with %d threads\n\n", NTHREADS);
    #endif

    fprintf(fp, "number of colors: %d\n", NCOLOR);
    fprintf(fp, "spacetime dimensionality: %d\n\n", STDIM);

    fprintf(fp, "lattice: %d", param->d_size[0]);
    for(i=1; i<STDIM; i++)
       {
       fprintf(fp, "x%d", param->d_size[i]);
       }
    fprintf(fp, "\n\n");

    fprintf(fp, "gfstep:        %lf\n", param->d_gfstep);
		fprintf(fp, "num_gfsteps    %d\n",  param->d_ngfsteps);
		fprintf(fp, "gf_meas_each   %d\n",  param->d_gf_meas_each);
    fprintf(fp, "\n");

    fprintf(fp, "randseed: %u\n", param->d_randseed);
    fprintf(fp, "\n");

    diff_sec = difftime(time_end, time_start);
    fprintf(fp, "Simulation time: %.3lf seconds\n", diff_sec );
    fprintf(fp, "\n");

    if(endian()==0)
      {
      fprintf(fp, "Little endian machine\n\n");
      }
    else
      {
      fprintf(fp, "Big endian machine\n\n");
      }

    fclose(fp);
    }

// print simulation parameters for the tracedef case
void print_parameters_tracedef(GParam const * const param, time_t time_start, time_t time_end, double acc)
    {
    FILE *fp;
    int i;
    double diff_sec;

    fp=fopen(param->d_log_file, "w");
    fprintf(fp, "+--------------------------------------------+\n");
    fprintf(fp, "| Simulation details for yang_mills_tracedef |\n");
    fprintf(fp, "+--------------------------------------------+\n\n");

    #ifdef OPENMP_MODE
     fprintf(fp, "using OpenMP with %d threads\n\n", NTHREADS);
    #endif

    fprintf(fp, "number of colors: %d\n", NCOLOR);
    fprintf(fp, "spacetime dimensionality: %d\n\n", STDIM);

    fprintf(fp, "lattice: %d", param->d_size[0]);
    for(i=1; i<STDIM; i++)
       {
       fprintf(fp, "x%d", param->d_size[i]);
       }
    fprintf(fp, "\n\n");

    fprintf(fp, "beta: %.10lf\n", param->d_beta);
    #ifdef THETA_MODE
      fprintf(fp, "theta: %.10lf\n", param->d_theta);
    #endif
    fprintf(fp, "h: %.10lf ", param->d_h[0]);
    for(i=1; i<(int) floor(NCOLOR/2.0); i++)
       {
       fprintf(fp, "%.10lf ", param->d_h[i]);
       }
    fprintf(fp, "\n\n");

    fprintf(fp, "sample:    %d\n", param->d_sample);
    fprintf(fp, "thermal:   %d\n", param->d_thermal);
    fprintf(fp, "overrelax: %d\n", param->d_overrelax);
    fprintf(fp, "measevery: %d\n", param->d_measevery);
    fprintf(fp, "\n");

    fprintf(fp, "start:                   %d\n", param->d_start);
    fprintf(fp, "saveconf_back_every:     %d\n", param->d_saveconf_back_every);
    fprintf(fp, "saveconf_analysis_every: %d\n", param->d_saveconf_analysis_every);
    fprintf(fp, "\n");

    fprintf(fp, "epsilon_metro: %.10lf\n", param->d_epsilon_metro);
    fprintf(fp, "metropolis acceptance: %.10lf\n", acc);
    fprintf(fp, "\n");

    fprintf(fp, "coolsteps:      %d\n", param->d_coolsteps);
    fprintf(fp, "coolrepeat:     %d\n", param->d_coolrepeat);
    fprintf(fp, "\n");

    fprintf(fp, "randseed: %u\n", param->d_randseed);
    fprintf(fp, "\n");

    diff_sec = difftime(time_end, time_start);
    fprintf(fp, "Simulation time: %.3lf seconds\n", diff_sec );
    fprintf(fp, "\n");

    if(endian()==0)
      {
      fprintf(fp, "Little endian machine\n\n");
      }
    else
      {
      fprintf(fp, "Big endian machine\n\n");
      }

    fclose(fp);
    }


// print simulation parameters
void print_parameters_tube_disc(GParam * param, time_t time_start, time_t time_end)
    {
    FILE *fp;
    int i;
    double diff_sec;

    fp=fopen(param->d_log_file, "w");
    fprintf(fp, "+---------------------------------------------+\n");
    fprintf(fp, "| Simulation details for yang_mills_tube_disc |\n");
    fprintf(fp, "+---------------------------------------------+\n\n");

    #ifdef OPENMP_MODE
     fprintf(fp, "using OpenMP with %d threads\n\n", NTHREADS);
    #endif

    fprintf(fp, "number of colors: %d\n", NCOLOR);
    fprintf(fp, "spacetime dimensionality: %d\n\n", STDIM);

    fprintf(fp, "lattice: %d", param->d_size[0]);
    for(i=1; i<STDIM; i++)
       {
       fprintf(fp, "x%d", param->d_size[i]);
       }
    fprintf(fp, "\n\n");

    fprintf(fp, "beta: %.10lf\n", param->d_beta);
    #ifdef THETA_MODE
      fprintf(fp, "theta: %.10lf\n", param->d_theta);
    #endif
    fprintf(fp, "\n");

    fprintf(fp, "sample:    %d\n", param->d_sample);
    fprintf(fp, "thermal:   %d\n", param->d_thermal);
    fprintf(fp, "overrelax: %d\n", param->d_overrelax);
    fprintf(fp, "measevery: %d\n", param->d_measevery);
    fprintf(fp, "\n");

    fprintf(fp, "start:                   %d\n", param->d_start);
    fprintf(fp, "saveconf_back_every:     %d\n", param->d_saveconf_back_every);
    fprintf(fp, "\n");

    fprintf(fp, "multihit:   %d\n", param->d_multihit);
    fprintf(fp, "levels for multileves: %d\n", NLEVELS);
    fprintf(fp, "multilevel steps: ");
    for(i=0; i<NLEVELS; i++)
       {
       fprintf(fp, "%d ", param->d_ml_step[i]);
       }
    fprintf(fp, "\n");
    fprintf(fp, "updates for levels: ");
    for(i=0; i<NLEVELS; i++)
       {
       fprintf(fp, "%d ", param->d_ml_upd[i]);
       }
    fprintf(fp, "\n");
    fprintf(fp, "dist_poly:   %d\n", param->d_dist_poly);
    fprintf(fp, "transv_dist: %d\n", param->d_trasv_dist);
    fprintf(fp, "plaq_dir: %d %d\n", param->d_plaq_dir[0], param->d_plaq_dir[1]);
    fprintf(fp, "\n");

    fprintf(fp, "randseed: %u\n", param->d_randseed);
    fprintf(fp, "\n");

    diff_sec = difftime(time_end, time_start);
    fprintf(fp, "Simulation time: %.3lf seconds\n", diff_sec );
    fprintf(fp, "\n");

    if(endian()==0)
      {
      fprintf(fp, "Little endian machine\n\n");
      }
    else
      {
      fprintf(fp, "Big endian machine\n\n");
      }

    fclose(fp);
    }


// print simulation parameters
void print_parameters_tube_conn(GParam * param, time_t time_start, time_t time_end)
    {
    FILE *fp;
    int i;
    double diff_sec;

    fp=fopen(param->d_log_file, "w");
    fprintf(fp, "+---------------------------------------------+\n");
    fprintf(fp, "| Simulation details for yang_mills_tube_conn |\n");
    fprintf(fp, "+---------------------------------------------+\n\n");

    #ifdef OPENMP_MODE
     fprintf(fp, "using OpenMP with %d threads\n\n", NTHREADS);
    #endif

    fprintf(fp, "number of colors: %d\n", NCOLOR);
    fprintf(fp, "spacetime dimensionality: %d\n\n", STDIM);

    fprintf(fp, "lattice: %d", param->d_size[0]);
    for(i=1; i<STDIM; i++)
       {
       fprintf(fp, "x%d", param->d_size[i]);
       }
    fprintf(fp, "\n\n");

    fprintf(fp, "beta: %.10lf\n", param->d_beta);
    #ifdef THETA_MODE
      fprintf(fp, "theta: %.10lf\n", param->d_theta);
    #endif
    fprintf(fp, "\n");

    fprintf(fp, "sample:    %d\n", param->d_sample);
    fprintf(fp, "thermal:   %d\n", param->d_thermal);
    fprintf(fp, "overrelax: %d\n", param->d_overrelax);
    fprintf(fp, "measevery: %d\n", param->d_measevery);
    fprintf(fp, "\n");

    fprintf(fp, "start:                   %d\n", param->d_start);
    fprintf(fp, "saveconf_back_every:     %d\n", param->d_saveconf_back_every);
    fprintf(fp, "\n");

    fprintf(fp, "multihit:   %d\n", param->d_multihit);
    fprintf(fp, "levels for multileves: %d\n", NLEVELS);
    fprintf(fp, "multilevel steps: ");
    for(i=0; i<NLEVELS; i++)
       {
       fprintf(fp, "%d ", param->d_ml_step[i]);
       }
    fprintf(fp, "\n");
    fprintf(fp, "updates for levels: ");
    for(i=0; i<NLEVELS; i++)
       {
       fprintf(fp, "%d ", param->d_ml_upd[i]);
       }
    fprintf(fp, "\n");
    fprintf(fp, "dist_poly:   %d\n", param->d_dist_poly);
    fprintf(fp, "transv_dist: %d\n", param->d_trasv_dist);
    fprintf(fp, "plaq_dir: %d %d\n", param->d_plaq_dir[0], param->d_plaq_dir[1]);
    fprintf(fp, "\n");

    fprintf(fp, "randseed: %u\n", param->d_randseed);
    fprintf(fp, "\n");

    diff_sec = difftime(time_end, time_start);
    fprintf(fp, "Simulation time: %.3lf seconds\n", diff_sec );
    fprintf(fp, "\n");

    if(endian()==0)
      {
      fprintf(fp, "Little endian machine\n\n");
      }
    else
      {
      fprintf(fp, "Big endian machine\n\n");
      }

    fclose(fp);
    }


// print simulation parameters
void print_parameters_tube_conn_long(GParam * param, time_t time_start, time_t time_end)
    {
    FILE *fp;
    int i;
    double diff_sec;

    fp=fopen(param->d_log_file, "w");
    fprintf(fp, "+--------------------------------------------------+\n");
    fprintf(fp, "| Simulation details for yang_mills_tube_conn_long |\n");
    fprintf(fp, "+--------------------------------------------------+\n\n");

    #ifdef OPENMP_MODE
     fprintf(fp, "using OpenMP with %d threads\n\n", NTHREADS);
    #endif

    fprintf(fp, "number of colors: %d\n", NCOLOR);
    fprintf(fp, "spacetime dimensionality: %d\n\n", STDIM);

    fprintf(fp, "lattice: %d", param->d_size[0]);
    for(i=1; i<STDIM; i++)
       {
       fprintf(fp, "x%d", param->d_size[i]);
       }
    fprintf(fp, "\n\n");

    fprintf(fp, "beta: %.10lf\n", param->d_beta);
    #ifdef THETA_MODE
      fprintf(fp, "theta: %.10lf\n", param->d_theta);
    #endif
    fprintf(fp, "\n");

    fprintf(fp, "sample:    %d\n", param->d_sample);
    fprintf(fp, "thermal:   %d\n", param->d_thermal);
    fprintf(fp, "overrelax: %d\n", param->d_overrelax);
    fprintf(fp, "measevery: %d\n", param->d_measevery);
    fprintf(fp, "\n");

    fprintf(fp, "start:                   %d\n", param->d_start);
    fprintf(fp, "saveconf_back_every:     %d\n", param->d_saveconf_back_every);
    fprintf(fp, "\n");

    fprintf(fp, "multihit:   %d\n", param->d_multihit);
    fprintf(fp, "levels for multileves: %d\n", NLEVELS);
    fprintf(fp, "multilevel steps: ");
    for(i=0; i<NLEVELS; i++)
       {
       fprintf(fp, "%d ", param->d_ml_step[i]);
       }
    fprintf(fp, "\n");
    fprintf(fp, "updates for levels: ");
    for(i=0; i<NLEVELS; i++)
       {
       fprintf(fp, "%d ", param->d_ml_upd[i]);
       }
    fprintf(fp, "\n");
    fprintf(fp, "level0_repeat:   %d\n", param->d_ml_level0_repeat);
    fprintf(fp, "dist_poly:   %d\n", param->d_dist_poly);
    fprintf(fp, "transv_dist: %d\n", param->d_trasv_dist);
    fprintf(fp, "plaq_dir: %d %d\n", param->d_plaq_dir[0], param->d_plaq_dir[1]);
    fprintf(fp, "\n");

    fprintf(fp, "randseed: %u\n", param->d_randseed);
    fprintf(fp, "\n");

    diff_sec = difftime(time_end, time_start);
    fprintf(fp, "Simulation time: %.3lf seconds\n", diff_sec );
    fprintf(fp, "\n");

    if(endian()==0)
      {
      fprintf(fp, "Little endian machine\n\n");
      }
    else
      {
      fprintf(fp, "Big endian machine\n\n");
      }

    fclose(fp);
    }



#endif

