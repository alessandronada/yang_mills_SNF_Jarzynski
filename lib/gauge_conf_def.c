#ifndef GAUGE_CONF_DEF_C
#define GAUGE_CONF_DEF_C

#include"../include/macro.h"

#ifdef HASH_MODE
  #include<openssl/md5.h>
#endif
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>

#include"../include/function_pointers.h"
#include"../include/gparam.h"
#include"../include/geometry.h"
#include"../include/gauge_conf.h"
#include"../include/tens_prod.h"

void init_gauge_conf_from_file_with_name(Gauge_Conf *GC, GParam const * const param, char const * const filename)
  {
  long r, j;
  int err;

  // allocate the local lattice
  err=posix_memalign((void**) &(GC->lattice), (size_t) DOUBLE_ALIGN, (size_t) param->d_volume * sizeof(GAUGE_GROUP *));
  if(err!=0)
    {
    fprintf(stderr, "Problems in allocating the lattice! (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  for(r=0; r<(param->d_volume); r++)
     {
     err=posix_memalign((void**) &(GC->lattice[r]), (size_t) DOUBLE_ALIGN, (size_t) STDIM * sizeof(GAUGE_GROUP));
     if(err!=0)
       {
       fprintf(stderr, "Problems in allocating the lattice! (%s, %d)\n", __FILE__, __LINE__);
       exit(EXIT_FAILURE);
       }
     }

  #ifdef THETA_MODE
    alloc_clover_array(GC, param);
  #endif

  // initialize lattice
  if(param->d_start==0) // ordered start
    {
    GAUGE_GROUP aux1, aux2;
    one(&aux1);

    GC->update_index=0;
    GC->evolution_index=0;

    for(r=0; r<(param->d_volume); r++)
       {
       for(j=0; j<STDIM; j++)
          {
          rand_matrix(&aux2);
          times_equal_real(&aux2, 0.001);
          plus_equal(&aux2, &aux1);
          unitarize(&aux2);
          equal_dag(&(GC->lattice[r][j]), &aux2);
          }
       }
    }
  if(param->d_start==1)  // random start
    {
    GAUGE_GROUP aux1;

    GC->update_index=0;
    GC->evolution_index=0;

    for(r=0; r<(param->d_volume); r++)
       {
       for(j=0; j<STDIM; j++)
          {
          rand_matrix(&aux1);
          equal(&(GC->lattice[r][j]), &aux1);
          }
       }
    }

  if(param->d_start==2) // initialize from stored conf
    {
    read_gauge_conf_from_file_with_name(GC, param, filename);
    }
  }

void init_gauge_conf(Gauge_Conf *GC, GParam const * const param)
	{
	init_gauge_conf_from_file_with_name(GC, param, param->d_conf_file);
	}

// used to allocate all replicas in the parallel tempering
void init_gauge_conf_replica(Gauge_Conf **GC, GParam  const * const param)
	{
	int i, err;
	
	// allocate the vector to store replicas
  err=posix_memalign((void **) GC, (size_t) DOUBLE_ALIGN, (size_t) param->d_N_replica_pt * sizeof(Gauge_Conf));
  if(err!=0)
    {
    fprintf(stderr, "Problems in allocating the parallel tempering replicas! (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }

  #ifdef OPENMP_MODE
  #pragma omp parallel for num_threads(NTHREADS) private(i)
  #endif
	for(i=0; i<param->d_N_replica_pt; i++)
		{
  	char filename[STD_STRING_LENGTH], replica_index[STD_STRING_LENGTH];
		strcpy(filename,param->d_conf_file); // filename = param->d_conf_file
		strcat(filename,"_replica_");			
		sprintf(replica_index, "%d", i);
		strcat(filename, replica_index); // filename = param->d_conf_file + "_replica_${i}"
		init_gauge_conf_from_file_with_name(&((*GC)[i]), param, filename);
		init_bound_cond(&((*GC)[i]),param,i);
		((*GC)[i]).conf_label=i;		
		}
	}

// initialization of the defect for the replica with label i
void init_bound_cond(Gauge_Conf *GC, GParam const * const param, int const i)
  {
	const int TRUE  = 1;
	const int FALSE = 0;
	long r;
	int err,j, is_on_defect;
	int cartcoord[STDIM];
	
	// for each value of defect_dir, determine the three orthogonal directions to it
	int perp_dir[4][3] = { {1, 2, 3}, {0, 2, 3}, {0, 1, 3}, {0, 1, 2} };
	
  //allocation of C[r][j]
  err=posix_memalign((void**) &(GC->C), (size_t) DOUBLE_ALIGN, (size_t) param->d_volume * sizeof(double *));
  if(err!=0)
    {
    fprintf(stderr, "Problems in allocating the defect! (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  for(r=0; r<(param->d_volume); r++)
     {
     err=posix_memalign((void**) &(GC->C[r]), (size_t) DOUBLE_ALIGN, (size_t) STDIM * sizeof(double));
     if(err!=0)
       {
       fprintf(stderr, "Problems in allocating the defect! (%s, %d)\n", __FILE__, __LINE__);
       exit(EXIT_FAILURE);
       }
     }

	// initialization of C[r][j]

	// start initializing them to 1
	for(r=0; r<param->d_volume; r++)
		for(j=0; j<STDIM; j++)
			{
			GC->C[r][j]=1.0;
			}

	// if there is more than 1 replica, initialize c(r) != 1 on the defect and along direction <defect_dir>
	if (param->d_N_replica_pt > 1)
	{
	j=param->d_defect_dir;
	for(r=0; r<param->d_volume; r++)
		{
		// check if r is on the defect or not
		si_to_cart(cartcoord, r, param);
		if( (cartcoord[param->d_defect_dir]==((param->d_size[param->d_defect_dir])-1)) &&
			(cartcoord[perp_dir[param->d_defect_dir][0]]<(param->d_L_defect[0])) &&
			(cartcoord[perp_dir[param->d_defect_dir][1]]<(param->d_L_defect[1])) &&
			(cartcoord[perp_dir[param->d_defect_dir][2]]<(param->d_L_defect[2])) ) is_on_defect=TRUE;
		else is_on_defect=FALSE;

		// if r is on defect assign bound cond
		if (is_on_defect==TRUE) {GC->C[r][j]=param->d_pt_bound_cond_coeff[i];}
		}
	}
  }

// initialize the defect of a configuration with a specific value of C
void init_single_conf_bound_cond(Gauge_Conf *GC, GParam const * const param, double const Cval)
{
	long r;
	int err;

	//allocation of C[r][j]
	err = posix_memalign((void**) &(GC->C), (size_t)DOUBLE_ALIGN, (size_t)param->d_volume * sizeof(double *));
	if (err != 0)
	{
		fprintf(stderr, "Problems in allocating the defect! (%s, %d)\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
	}
	for (r = 0; r < (param->d_volume); r++)
	{
		err = posix_memalign((void**) &(GC->C[r]), (size_t)DOUBLE_ALIGN, (size_t)STDIM * sizeof(double));
		if (err != 0)
		{
			fprintf(stderr, "Problems in allocating the defect! (%s, %d)\n", __FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}
	}

	set_bound_cond(GC, param, Cval);
}

// set the defect of a configuration with a specific value of C
void set_bound_cond(Gauge_Conf *GC, GParam const * const param, double const Cval)
{
	const int TRUE = 1;
	const int FALSE = 0;
	long r;
	int j, is_on_defect;
	int cartcoord[STDIM];

	// for each value of defect_dir, determine the three orthogonal directions to it
	int perp_dir[4][3] = { {1, 2, 3}, {0, 2, 3}, {0, 1, 3}, {0, 1, 2} };

	// initialization of C[r][j]

	// start initializing them to 1
	for (r = 0; r < param->d_volume; r++)
		for (j = 0; j < STDIM; j++)
		{
			GC->C[r][j] = 1.0;
		}

	j = param->d_defect_dir;
	for (r = 0; r < param->d_volume; r++)
	{
		// check if r is on the defect or not
		si_to_cart(cartcoord, r, param);
		if ((cartcoord[param->d_defect_dir] == ((param->d_size[param->d_defect_dir]) - 1)) &&
			(cartcoord[perp_dir[param->d_defect_dir][0]] < (param->d_L_defect[0])) &&
			(cartcoord[perp_dir[param->d_defect_dir][1]] < (param->d_L_defect[1])) &&
			(cartcoord[perp_dir[param->d_defect_dir][2]] < (param->d_L_defect[2]))) is_on_defect = TRUE;
		else is_on_defect = FALSE;

		// if r is on defect assign bound cond
		if (is_on_defect == TRUE) { GC->C[r][j] = Cval; }
	}
}

void read_gauge_conf_from_file_with_name(Gauge_Conf *GC, GParam const * const param, char const * const filename)
  {
  FILE *fp;
  int i, dimension, tmp_i;
  int err, mu;
  long lex, si;
  GAUGE_GROUP matrix;
  #ifdef HASH_MODE
    char md5sum_new[2*MD5_DIGEST_LENGTH+1];
    char md5sum_old[2*MD5_DIGEST_LENGTH+1];
  #else
    char md5sum_old[2*STD_STRING_LENGTH+1]={0};
  #endif

  fp=fopen(filename, "r"); // open the configuration file
  if(fp==NULL)
    {
    fprintf(stderr, "Error in opening the file %s (%s, %d)\n", filename, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  else // read the txt header of the configuration
    {
    err=fscanf(fp, "%d", &dimension);
    if(err!=1)
      {
      fprintf(stderr, "Error in reading the file %s (%s, %d)\n", filename, __FILE__, __LINE__);
      exit(EXIT_FAILURE);
      }
    if(dimension != STDIM)
      {
      fprintf(stderr, "The space time dimension of the configuration (%d) does not coincide with the one of the global parameter (%d)\n",
              dimension, STDIM);
      exit(EXIT_FAILURE);
      }

    for(i=0; i<STDIM; i++)
       {
       err=fscanf(fp, "%d", &tmp_i);
       if(err!=1)
         {
         fprintf(stderr, "Error in reading the file %s (%s, %d)\n", filename, __FILE__, __LINE__);
         exit(EXIT_FAILURE);
         }
       if(tmp_i != param->d_size[i])
         {
         fprintf(stderr, "The size of the configuration lattice does not coincide with the one of the global parameter\n");
         exit(EXIT_FAILURE);
         }
       }

    err=fscanf(fp, "%ld %ld %s\n", &(GC->update_index), &(GC->evolution_index), md5sum_old);
    if(err!=3)
      {
      fprintf(stderr, "Error in reading the file %s (%s, %d)\n", filename, __FILE__, __LINE__);
      exit(EXIT_FAILURE);
      }

    fclose(fp);
    }

  fp=fopen(filename, "rb"); // open the configuration file in binary
  if(fp==NULL)
    {
    fprintf(stderr, "Error in opening the file %s (%s, %d)\n", filename, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  else
    {
    // read again the header
    err=0;
    while(err!='\n')
         {
         err=fgetc(fp);
         }

    for(lex=0; lex<param->d_volume; lex++)
       {
       si=lex_to_si(lex, param);
       for(mu=0; mu<STDIM; mu++)
          {
          read_from_binary_file_bigen(fp, &matrix);

          equal(&(GC->lattice[si][mu]), &matrix);
          }
       }
    fclose(fp);

    #ifdef HASH_MODE
      // compute the new md5sum and check for consistency
      compute_md5sum_conf(md5sum_new, GC, param);
      if(strncmp(md5sum_old, md5sum_new, 2*MD5_DIGEST_LENGTH+1)!=0)
        {
        fprintf(stderr, "The computed md5sum %s does not match the stored %s (%s, %d)\n", md5sum_new, md5sum_old, __FILE__, __LINE__);
        exit(EXIT_FAILURE);
        }
    #endif
    }
  }

void free_gauge_conf(Gauge_Conf *GC, GParam const * const param)
  {
  long i;

  for(i=0; i<(param->d_volume); i++)
     {
     free(GC->lattice[i]);
     }
  free(GC->lattice);

  #ifdef THETA_MODE
    end_clover_array(GC, param);
  #endif
  }
	
void free_replica(Gauge_Conf *GC, GParam const * const param)
	{
	int i;
	for(i=0; i<param->d_N_replica_pt; i++)
		{
		free_gauge_conf(&(GC[i]), param);
		free_bound_cond(&(GC[i]), param);
		}
	free(GC);
	}
	
void free_bound_cond(Gauge_Conf *GC, GParam const * const param)
  {
  long r;
  for(r=0; r<(param->d_volume); r++)
     {
     free(GC->C[r]);
     }
  free(GC->C);
  }


// save a configuration in ILDG-like format
void write_conf_on_file_with_name(Gauge_Conf const * const GC,
                            GParam const * const param,
                            char const * const namefile)
  {
  long si, lex;
  int i, mu;
  #ifdef HASH_MODE
    char md5sum[2*MD5_DIGEST_LENGTH+1];
  #else
     char md5sum[2*STD_STRING_LENGTH+1]={0};
  #endif
  FILE *fp;

  #ifdef HASH_MODE
    compute_md5sum_conf(md5sum, GC, param);
  #endif

  fp=fopen(namefile, "w"); // open the configuration file
  if(fp==NULL)
    {
    fprintf(stderr, "Error in opening the file %s (%s, %d)\n", namefile, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  else
    {
    fprintf(fp, "%d ", STDIM);
    for(i=0; i<STDIM; i++)
       {
       fprintf(fp, "%d ", param->d_size[i]);
       }
    fprintf(fp, "%ld %ld %s\n", GC->update_index, GC->evolution_index, md5sum);
    }
  fclose(fp);

  fp=fopen(namefile, "ab"); // open the configuration file in binary mode
  if(fp==NULL)
    {
    fprintf(stderr, "Error in opening the file %s (%s, %d)\n", namefile, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  else
    {
    for(lex=0; lex<param->d_volume; lex++)
       {
       si=lex_to_si(lex, param);
       for(mu=0; mu<STDIM; mu++)
          {
          print_on_binary_file_bigen(fp, &(GC->lattice[si][mu]) );
          }
       }
    fclose(fp);
    }
  }


void write_conf_on_file(Gauge_Conf const * const GC, GParam const * const param)
  {
  write_conf_on_file_with_name(GC, param, param->d_conf_file);
  }

void write_evolution_conf_on_file(Gauge_Conf const* const GC, GParam const* const param, int obc)
{	
    char filename[STD_STRING_LENGTH], evolution_index[STD_STRING_LENGTH];
    strcpy(filename, param->d_conf_file);
    if (obc != 0)
        strcat(filename, "_ev_");
    else
        strcat(filename, "_PBC_ev_");
    sprintf(evolution_index, "%ld", GC->evolution_index);
    strcat(filename, evolution_index); // filename = d_conf_file + "_ev_${i}"
    write_conf_on_file_with_name(GC, param, filename);
}

void write_replica_on_file(Gauge_Conf const * const GC, GParam const * const param)
	{
	int i;

	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS) private(i)
	#endif	
	for(i=0; i<param->d_N_replica_pt; i++)
	{
		char filename[STD_STRING_LENGTH], replica_index[STD_STRING_LENGTH];
		strcpy(filename,param->d_conf_file);
		strcat(filename,"_replica_");
		sprintf(replica_index, "%d", i);
		strcat(filename,replica_index); // filename = d_conf_file + "_replica_${i}"
		write_conf_on_file_with_name(&(GC[i]),param,filename);
		}
	}
	
void write_replica_on_file_back(Gauge_Conf const * const GC, GParam const * const param)
	{
	int i;
	static int counter=0;
	
	#ifdef OPENMP_MODE
	#pragma omp parallel for num_threads(NTHREADS) private(i)
	#endif
	for(i=0; i<param->d_N_replica_pt; i++)
	{
		char filename[STD_STRING_LENGTH], replica_index[STD_STRING_LENGTH], aux_back[STD_STRING_LENGTH];
		if(counter==0)
			sprintf(aux_back, "_back0");
		else
			sprintf(aux_back, "_back1");
		strcpy(filename,param->d_conf_file);
		strcat(filename,"_replica_");
		sprintf(replica_index, "%d", i);
		strcat(filename, replica_index);
		strcat(filename, aux_back); // filename = d_conf_file + "_replica_${i}_back${counter}"
		write_conf_on_file_with_name(&(GC[i]),param,filename);
	}
	counter=1-counter;
	}

void write_conf_on_file_back(Gauge_Conf const * const GC, GParam const * const param)
  {
  char name[STD_STRING_LENGTH], aux[STD_STRING_LENGTH];
  static int counter=0;

  strcpy(name, param->d_conf_file);
  if(counter==0)
    {
    sprintf(aux, "_back0");
    }
  else
    {
    sprintf(aux, "_back1");
    }
  strcat(name, aux);

  write_conf_on_file_with_name(GC, param, name);

  counter=1-counter;
  }


// allocate GC and initialize with GC2
void init_gauge_conf_from_gauge_conf(Gauge_Conf *GC, Gauge_Conf const * const GC2, GParam const * const param) 
  {
  long r;
  int mu, err;

  // allocate the lattice
  err=posix_memalign((void**)&(GC->lattice), (size_t) DOUBLE_ALIGN, (size_t) param->d_volume * sizeof(GAUGE_GROUP *));
  if(err!=0)
    {
    fprintf(stderr, "Problems in allocating the lattice! (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  for(r=0; r<(param->d_volume); r++)
     {
     err=posix_memalign((void**)&(GC->lattice[r]), (size_t) DOUBLE_ALIGN, (size_t) STDIM * sizeof(GAUGE_GROUP));
     if(err!=0)
       {
       fprintf(stderr, "Problems in allocating the lattice! (%s, %d)\n", __FILE__, __LINE__);
       exit(EXIT_FAILURE);
       }
     }

  #ifdef THETA_MODE
    alloc_clover_array(GC, param);
  #endif

  // initialize GC
  for(r=0; r<(param->d_volume); r++)
     {
     for(mu=0; mu<STDIM; mu++)
        {
        equal(&(GC->lattice[r][mu]), &(GC2->lattice[r][mu]) );
        }
     }

  GC->update_index=GC2->update_index;
  GC->evolution_index=GC2->evolution_index;
  }

// allocate GC and initialize with GC2
void copy_gauge_conf_from_gauge_conf(Gauge_Conf *GC, Gauge_Conf const * const GC2, GParam const * const param)
{
	long r;
	int mu; //, err;

	// initialize GC
	for (r = 0; r < (param->d_volume); r++)
	{
		for (mu = 0; mu < STDIM; mu++)
		{
			equal(&(GC->lattice[r][mu]), &(GC2->lattice[r][mu]));
		}
	}

	GC->update_index = GC2->update_index;
	GC->evolution_index=GC2->evolution_index;
}


// compute the md5sum of the configuration and save it in res, that is a char[2*MD5_DIGEST_LENGTH]
void compute_md5sum_conf(char *res, Gauge_Conf const * const GC, GParam const * const param)
  {
  #ifdef HASH_MODE
    MD5_CTX mdContext;
    unsigned char c[MD5_DIGEST_LENGTH];
    long si, lex;
    GAUGE_GROUP matrix;
    int mu, k;

    MD5_Init(&mdContext);
    for(lex=0; lex<param->d_volume; lex++)
       {
       si=lex_to_si(lex, param);
       for(mu=0; mu<STDIM; mu++)
          {
          equal(&matrix, &(GC->lattice[si][mu]));

          #if NCOLOR==1
            MD5_Update(&mdContext, &(matrix.comp), sizeof(double complex));
          #elif NCOLOR==2
            for(k=0; k<4; k++)
               {
               MD5_Update(&mdContext, &(matrix.comp[k]), sizeof(double));
               }
          #else
            for(k=0; k<NCOLOR*NCOLOR; k++)
               {
               MD5_Update(&mdContext, &(matrix.comp[k]), sizeof(double complex));
               }
          #endif
          }
       }
    MD5_Final(c, &mdContext);

    for(k = 0; k < MD5_DIGEST_LENGTH; k++)
       {
       sprintf(&(res[2*k]), "%02x", c[k]);
       }
  #else
    // just to avoid warning at compile time
    (void) res;
    (void) GC;
    (void) param;
  #endif
  }


// allocate the ml_polycorr arrays and related stuff
void alloc_polycorr_stuff(Gauge_Conf *GC,
                          GParam const * const param)
  {
  int i, j, err;

  err=posix_memalign((void**)&(GC->ml_polycorr), (size_t) DOUBLE_ALIGN, (size_t) NLEVELS * sizeof(TensProd **));
  if(err!=0)
    {
    fprintf(stderr, "Problems in allocating ml_polycorr (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  else
    {
    for(i=0; i<NLEVELS; i++)
       {
       err=posix_memalign((void**)&(GC->ml_polycorr[i]), (size_t) DOUBLE_ALIGN, (size_t) (param->d_size[0] / param->d_ml_step[i]) * sizeof(TensProd *));
       if(err!=0)
         {
         fprintf(stderr, "Problems in allocating ml_polycorr[%d] (%s, %d)\n", i, __FILE__, __LINE__);
         exit(EXIT_FAILURE);
         }
       else
         {
         for(j=0; j<(param->d_size[0]/param->d_ml_step[i]); j++)
            {
            err=posix_memalign((void**)&(GC->ml_polycorr[i][j]), (size_t) DOUBLE_ALIGN, (size_t) param->d_space_vol * sizeof(TensProd));
            if(err!=0)
              {
              fprintf(stderr, "Problems in allocating ml_polycorr[%d][%d] (%s, %d)\n", i, j, __FILE__, __LINE__);
              exit(EXIT_FAILURE);
              }
            }
         }
       }
    }

  err=posix_memalign((void**)&(GC->loc_poly), (size_t) DOUBLE_ALIGN, (size_t) (param->d_size[0]/param->d_ml_step[NLEVELS-1]) * sizeof(GAUGE_GROUP *));
  if(err!=0)
    {
    fprintf(stderr, "Problems in allocating loc_poly (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  else
    {
    for(i=0; i<param->d_size[0]/param->d_ml_step[NLEVELS-1]; i++)
       {
       err=posix_memalign((void**)&(GC->loc_poly[i]), (size_t) DOUBLE_ALIGN, (size_t) param->d_space_vol *sizeof(GAUGE_GROUP));
       if(err!=0)
         {
         fprintf(stderr, "Problems in allocating loc_poly (%s, %d)\n", __FILE__, __LINE__);
         exit(EXIT_FAILURE);
         }
       }
    }
  }


// free the ml_polycorr arrays and related stuff
void free_polycorr_stuff(Gauge_Conf *GC,
                         GParam const * const param)
  {
  int i, j;

  for(i=0; i<NLEVELS; i++)
     {
     for(j=0; j<(param->d_size[0]/param->d_ml_step[i]); j++)
        {
        free(GC->ml_polycorr[i][j]);
        }
     free(GC->ml_polycorr[i]);
     }
  free(GC->ml_polycorr);

  for(i=0; i<param->d_size[0]/param->d_ml_step[NLEVELS-1]; i++)
     {
     free(GC->loc_poly[i]);
     }
  free(GC->loc_poly);
  }


// save ml_polycorr[0] arrays on file
void write_polycorr_on_file(Gauge_Conf const * const GC,
                            GParam const * const param,
                            int iteration)
  {
  long i;
  int j;
  #ifdef HASH_MODE
    char md5sum[2*MD5_DIGEST_LENGTH+1];
  #else
    char md5sum[2*STD_STRING_LENGTH+1]={0};
  #endif
  FILE *fp;

  #ifdef HASH_MODE
    compute_md5sum_polycorr(md5sum, GC, param);
  #endif

  fp=fopen(param->d_ml_file, "w"); // open the configuration file
  if(fp==NULL)
    {
    fprintf(stderr, "Error in opening the file %s (%s, %d)\n", param->d_ml_file, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  else
    {
    fprintf(fp, "%ld %d %s\n", param->d_space_vol, iteration, md5sum);
    }
  fclose(fp);

  fp=fopen(param->d_ml_file, "ab"); // open the configuration file in binary mode
  if(fp==NULL)
    {
    fprintf(stderr, "Error in opening the file %s (%s, %d)\n", param->d_ml_file, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  else
    {
    for(j=0; j<param->d_size[0]/param->d_ml_step[0]; j++)
       {
       for(i=0; i<(param->d_space_vol); i++)
          {
          print_on_binary_file_bigen_TensProd(fp, &(GC->ml_polycorr[0][j][i]));
          }
       }

    fclose(fp);
    }
  }


// read ml_polycorr[0] arrays from file
void read_polycorr_from_file(Gauge_Conf const * const GC,
                             GParam const * const param,
                             int *iteration)
  {
  long i, loc_space_vol;
  int j;
  FILE *fp;
  #ifdef HASH_MODE
    char md5sum_new[2*MD5_DIGEST_LENGTH+1];
    char md5sum_old[2*MD5_DIGEST_LENGTH+1];
  #else
    char md5sum_old[2*STD_STRING_LENGTH+1]={0};
  #endif

  fp=fopen(param->d_ml_file, "r"); // open the multilevel file
  if(fp==NULL)
    {
    fprintf(stderr, "Error in opening the file %s (%s, %d)\n", param->d_ml_file, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  else
    {
    i=fscanf(fp, "%ld %d %s\n", &loc_space_vol, iteration, md5sum_old);
    if(i!=3)
      {
      fprintf(stderr, "Error in reading the file %s (%s, %d)\n", param->d_ml_file, __FILE__, __LINE__);
      exit(EXIT_FAILURE);
      }
    if(loc_space_vol != param->d_space_vol)
      {
      fprintf(stderr, "Error: space_vol in the multilevel file %s is different from the one in the input (%s, %d)\n",
              param->d_ml_file, __FILE__, __LINE__);
      exit(EXIT_FAILURE);
      }
    }
  fclose(fp);

  fp=fopen(param->d_ml_file, "rb"); // open the multilevel file in binary mode
  if(fp==NULL)
    {
    fprintf(stderr, "Error in opening the file %s (%s, %d)\n", param->d_ml_file, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  else
    {
    // read again the header: loc_space_vol, iteration, hash
    i=0;
    while(i!='\n')
         {
         i=fgetc(fp);
         }

    for(j=0; j<param->d_size[0]/param->d_ml_step[0]; j++)
       {
       for(i=0; i<(param->d_space_vol); i++)
          {
          read_from_binary_file_bigen_TensProd(fp, &(GC->ml_polycorr[0][j][i]));
          }
       }

    fclose(fp);
    }

  #ifdef HASH_MODE
    // compute the new md5sum and check for consistency
    compute_md5sum_polycorr(md5sum_new, GC, param);
    if(strncmp(md5sum_old, md5sum_new, 2*MD5_DIGEST_LENGTH+1)!=0)
      {
      fprintf(stderr, "The computed md5sum %s of the multilevel file does not match the stored %s\n", md5sum_new, md5sum_old);
      exit(EXIT_FAILURE);
      }
  #endif
  }


// compute the md5sum of the ml_polycorr[0] arrays and save it in res, that is a char[2*MD5_DIGEST_LENGTH]
void compute_md5sum_polycorr(char *res, Gauge_Conf const * const GC, GParam const * const param)
  {
  #ifdef HASH_MODE
    MD5_CTX mdContext;
    unsigned char c[MD5_DIGEST_LENGTH];
    long i;
    int j;
    int n1, n2, n3, n4;

    MD5_Init(&mdContext);

    for(j=0; j<param->d_size[0]/param->d_ml_step[0]; j++)
       {
       for(i=0; i<(param->d_space_vol); i++)
          {
          for(n1=0; n1<NCOLOR; n1++)
             {
             for(n2=0; n2<NCOLOR; n2++)
                {
                for(n3=0; n3<NCOLOR; n3++)
                   {
                   for(n4=0; n4<NCOLOR; n4++)
                      {
                      MD5_Update(&mdContext, &((GC->ml_polycorr[0][j][i]).comp[n1][n2][n3][n4]), sizeof(double complex));
                      }
                   }
                }
             }
          }
       }

    MD5_Final(c, &mdContext);

    for(i = 0; i < MD5_DIGEST_LENGTH; i++)
       {
       sprintf(&(res[2*i]), "%02x", c[i]);
       }
  #else
    // just to avoid warning at compile time
    (void) res;
    (void) GC;
    (void) param;
  #endif
  }


// allocate the ml_polycorradj arrays
void alloc_polycorradj(Gauge_Conf *GC,
                       GParam const * const param)
  {
  int i, j, err;

  err=posix_memalign((void**)&(GC->ml_polycorradj), (size_t) DOUBLE_ALIGN, (size_t) NLEVELS * sizeof(TensProdAdj **));
  if(err!=0)
    {
    fprintf(stderr, "Problems in allocating ml_polycorradj (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  else
    {
    for(i=0; i<NLEVELS; i++)
       {
       err=posix_memalign((void**)&(GC->ml_polycorradj[i]), (size_t) DOUBLE_ALIGN, (size_t) (param->d_size[0] / param->d_ml_step[i]) * sizeof(TensProdAdj *));
       if(err!=0)
         {
         fprintf(stderr, "Problems in allocating ml_polycorradj[%d] (%s, %d)\n", i, __FILE__, __LINE__);
         exit(EXIT_FAILURE);
         }
       else
         {
         for(j=0; j<(param->d_size[0]/param->d_ml_step[i]); j++)
            {
            err=posix_memalign((void**)&(GC->ml_polycorradj[i][j]), (size_t) DOUBLE_ALIGN, (size_t) param->d_space_vol * sizeof(TensProdAdj));
            if(err!=0)
              {
              fprintf(stderr, "Problems in allocating ml_polycorrag[%d][%d] (%s, %d)\n", i, j, __FILE__, __LINE__);
              exit(EXIT_FAILURE);
              }
            }
         }
       }
    }
  }


// free the ml_polycorradj arrays
void free_polycorradj(Gauge_Conf *GC,
                      GParam const * const param)
  {
  int i, j;

  for(i=0; i<NLEVELS; i++)
     {
     for(j=0; j<(param->d_size[0]/param->d_ml_step[i]); j++)
        {
        free(GC->ml_polycorradj[i][j]);
        }
     free(GC->ml_polycorradj[i]);
     }
  free(GC->ml_polycorradj);
  }








// allocate the ml_polycorr, polyplaq arrays and related stuff
void alloc_tube_disc_stuff(Gauge_Conf *GC,
                           GParam const * const param)
  {
  int i, err;

  alloc_polycorr_stuff(GC, param);

  err=posix_memalign((void**)&(GC->ml_polyplaq), (size_t) DOUBLE_ALIGN, (size_t) NLEVELS * sizeof(TensProd *));
  if(err!=0)
    {
    fprintf(stderr, "Problems in allocating ml_polyplaq (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  else
    {
    for(i=0; i<NLEVELS; i++)
       {
       err=posix_memalign((void**)&(GC->ml_polyplaq[i]), (size_t) DOUBLE_ALIGN, (size_t) param->d_space_vol * sizeof(TensProd));
       if(err!=0)
         {
         fprintf(stderr, "Problems in allocating ml_polyplaq[%d] (%s, %d)\n", i, __FILE__, __LINE__);
         exit(EXIT_FAILURE);
         }
       }
    }

  err=posix_memalign((void**)&(GC->loc_plaq), (size_t)DOUBLE_ALIGN, (size_t) param->d_space_vol *sizeof(double complex));
  if(err!=0)
    {
    fprintf(stderr, "Problems in allocating loc_plaq (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  }


// free the ml_polycorr, ml_polyplaq arrays and relates stuff
void free_tube_disc_stuff(Gauge_Conf *GC,
                          GParam const * const param)
  {
  int i;

  free_polycorr_stuff(GC, param);

  for(i=0; i<NLEVELS; i++)
     {
     free(GC->ml_polyplaq[i]);
     }
  free(GC->ml_polyplaq);

  free(GC->loc_plaq);
  }


// allocate the polycorr, polyplaq, polyplaqconn arrays and related stuff
void alloc_tube_conn_stuff(Gauge_Conf *GC,
                           GParam const * const param)
  {
  int i, err;

  alloc_tube_disc_stuff(GC, param);

  err=posix_memalign((void**)&(GC->ml_polyplaqconn), (size_t) DOUBLE_ALIGN, (size_t) NLEVELS * sizeof(TensProd *));
  if(err!=0)
    {
    fprintf(stderr, "Problems in allocating ml_polyplaqconn (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  else
    {
    for(i=0; i<NLEVELS; i++)
       {
       err=posix_memalign((void**)&(GC->ml_polyplaqconn[i]), (size_t) DOUBLE_ALIGN, (size_t) param->d_space_vol * sizeof(TensProd));
       if(err!=0)
         {
         fprintf(stderr, "Problems in allocating ml_polyplaqconn[%d] (%s, %d)\n", i, __FILE__, __LINE__);
         exit(EXIT_FAILURE);
         }
       }
    }

  err=posix_memalign((void**)&(GC->loc_plaqconn), (size_t)DOUBLE_ALIGN, (size_t) param->d_space_vol *sizeof(GAUGE_GROUP));
  if(err!=0)
    {
    fprintf(stderr, "Problems in allocating loc_polyplaqconn (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }

  }


// free the polycorr, polyplaq, polyplaqconn arrays and related stuff
void free_tube_conn_stuff(Gauge_Conf *GC,
                          GParam const * const param)
  {
  int i, j;

  for(i=0; i<NLEVELS; i++)
     {
     for(j=0; j<(param->d_size[0]/param->d_ml_step[i]); j++)
        {
        free(GC->ml_polycorr[i][j]);
        }
     free(GC->ml_polycorr[i]);
     }
  free(GC->ml_polycorr);

  free(GC->loc_plaqconn);
  }


// save ml_polycorr[0], ml_polyplaq[0] and ml_polyplaqconn[0] arrays on file
void write_tube_conn_stuff_on_file(Gauge_Conf const * const GC,
                                   GParam const * const param,
                                   int iteration)
  {
  long i;
  int j;
  #ifdef HASH_MODE
    char md5sum[2*MD5_DIGEST_LENGTH+1];
  #else
    char md5sum[2*STD_STRING_LENGTH+1]={0};
  #endif
  FILE *fp;

  #ifdef HASH_MODE
    compute_md5sum_tube_conn_stuff(md5sum, GC, param);
  #endif

  fp=fopen(param->d_ml_file, "w"); // open the configuration file
  if(fp==NULL)
    {
    fprintf(stderr, "Error in opening the file %s (%s, %d)\n", param->d_ml_file, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  else
    {
    fprintf(fp, "%ld %d %s\n", param->d_space_vol, iteration, md5sum);
    }
  fclose(fp);

  fp=fopen(param->d_ml_file, "ab"); // open the configuration file in binary mode
  if(fp==NULL)
    {
    fprintf(stderr, "Error in opening the file %s (%s, %d)\n", param->d_ml_file, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  else
    {
    for(j=0; j<param->d_size[0]/param->d_ml_step[0]; j++)
       {
       for(i=0; i<(param->d_space_vol); i++)
          {
          print_on_binary_file_bigen_TensProd(fp, &(GC->ml_polycorr[0][j][i]));
          }
       }
    for(i=0; i<(param->d_space_vol); i++)
       {
       print_on_binary_file_bigen_TensProd(fp, &(GC->ml_polyplaq[0][i]));
       }
    for(i=0; i<(param->d_space_vol); i++)
       {
       print_on_binary_file_bigen_TensProd(fp, &(GC->ml_polyplaqconn[0][i]));
       }

    fclose(fp);
    }
  }


// read ml_polycorr[0], ml_polyplaq[0] and ml_polyplaqconn[0] arrays from file
void read_tube_conn_stuff_from_file(Gauge_Conf const * const GC,
                                    GParam const * const param,
                                    int *iteration)
  {
  long i, loc_space_vol;
  int j;
  FILE *fp;
  #ifdef HASH_MODE
    char md5sum_new[2*MD5_DIGEST_LENGTH+1];
    char md5sum_old[2*MD5_DIGEST_LENGTH+1];
  #else
    char md5sum_old[2*STD_STRING_LENGTH+1]={0};
  #endif

  fp=fopen(param->d_ml_file, "r"); // open the multilevel file
  if(fp==NULL)
    {
    fprintf(stderr, "Error in opening the file %s (%s, %d)\n", param->d_ml_file, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  else
    {
    i=fscanf(fp, "%ld %d %s\n", &loc_space_vol, iteration, md5sum_old);
    if(i!=3)
      {
      fprintf(stderr, "Error in reading the file %s (%s, %d)\n", param->d_ml_file, __FILE__, __LINE__);
      exit(EXIT_FAILURE);
      }
    if(loc_space_vol != param->d_space_vol)
      {
      fprintf(stderr, "Error: space_vol in the multilevel file %s is different from the one in the input (%s, %d)\n",
              param->d_ml_file, __FILE__, __LINE__);
      exit(EXIT_FAILURE);
      }
    }
  fclose(fp);

  fp=fopen(param->d_ml_file, "rb"); // open the multilevel file in binary mode
  if(fp==NULL)
    {
    fprintf(stderr, "Error in opening the file %s (%s, %d)\n", param->d_ml_file, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  else
    {
    // read again the header: loc_space_vol, iteration, hash
    i=0;
    while(i!='\n')
         {
         i=fgetc(fp);
         }

    for(j=0; j<param->d_size[0]/param->d_ml_step[0]; j++)
       {
       for(i=0; i<(param->d_space_vol); i++)
          {
          read_from_binary_file_bigen_TensProd(fp, &(GC->ml_polycorr[0][j][i]));
          }
       }
    for(i=0; i<(param->d_space_vol); i++)
       {
       read_from_binary_file_bigen_TensProd(fp, &(GC->ml_polyplaq[0][i]));
       }
    for(i=0; i<(param->d_space_vol); i++)
       {
       read_from_binary_file_bigen_TensProd(fp, &(GC->ml_polyplaqconn[0][i]));
       }

    fclose(fp);
    }

  #ifdef HASH_MODE
    // compute the new md5sum and check for consistency
    compute_md5sum_tube_conn_stuff(md5sum_new, GC, param);
    if(strncmp(md5sum_old, md5sum_new, 2*MD5_DIGEST_LENGTH+1)!=0)
      {
      fprintf(stderr, "The computed md5sum %s of the multilevel file does not match the stored %s\n", md5sum_new, md5sum_old);
      exit(EXIT_FAILURE);
      }
  #endif
  }


// compute the md5sum of the ml_polycorr[0], ml_polyplaq[0] and ml_polyplaqconn[0] arrays and save it in res, that is a char[2*MD5_DIGEST_LENGTH]
void compute_md5sum_tube_conn_stuff(char *res, Gauge_Conf const * const GC, GParam const * const param)
  {
  #ifdef HASH_MODE
    MD5_CTX mdContext;
    unsigned char c[MD5_DIGEST_LENGTH];
    long i;
    int j;
    int n1, n2, n3, n4;

    MD5_Init(&mdContext);

    for(j=0; j<param->d_size[0]/param->d_ml_step[0]; j++)
       {
       for(i=0; i<(param->d_space_vol); i++)
          {
          for(n1=0; n1<NCOLOR; n1++)
             {
             for(n2=0; n2<NCOLOR; n2++)
                {
                for(n3=0; n3<NCOLOR; n3++)
                   {
                   for(n4=0; n4<NCOLOR; n4++)
                      {
                      MD5_Update(&mdContext, &((GC->ml_polycorr[0][j][i]).comp[n1][n2][n3][n4]), sizeof(double complex));
                      }
                   }
                }
             }
          }
       }

    for(i=0; i<(param->d_space_vol); i++)
       {
       for(n1=0; n1<NCOLOR; n1++)
          {
          for(n2=0; n2<NCOLOR; n2++)
             {
             for(n3=0; n3<NCOLOR; n3++)
                {
                for(n4=0; n4<NCOLOR; n4++)
                   {
                   MD5_Update(&mdContext, &((GC->ml_polyplaq[0][i]).comp[n1][n2][n3][n4]), sizeof(double complex));
                   }
                }
             }
          }
       }

    for(i=0; i<(param->d_space_vol); i++)
       {
       for(n1=0; n1<NCOLOR; n1++)
          {
          for(n2=0; n2<NCOLOR; n2++)
             {
             for(n3=0; n3<NCOLOR; n3++)
                {
                for(n4=0; n4<NCOLOR; n4++)
                   {
                   MD5_Update(&mdContext, &((GC->ml_polyplaqconn[0][i]).comp[n1][n2][n3][n4]), sizeof(double complex));
                   }
                }
             }
          }
       }

    MD5_Final(c, &mdContext);

    for(i = 0; i < MD5_DIGEST_LENGTH; i++)
       {
       sprintf(&(res[2*i]), "%02x", c[i]);
       }
  #else
    // just to avoid warning at compile time
    (void) res;
    (void) GC;
    (void) param;
  #endif
  }


// allocate the clovers arrays
void alloc_clover_array(Gauge_Conf *GC,
                       GParam const * const param)
  {
  int i, err;
  long r;

  err=posix_memalign((void**)&(GC->clover_array), (size_t)DOUBLE_ALIGN, (size_t) param->d_volume *sizeof(GAUGE_GROUP **));
  if(err!=0)
    {
    fprintf(stderr, "Problems in allocating clovers (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  else
    {
    for(r=0; r<param->d_volume; r++)
       {
       err=posix_memalign((void**)&(GC->clover_array[r]), (size_t)DOUBLE_ALIGN, STDIM*sizeof(GAUGE_GROUP *));
       if(err!=0)
         {
         fprintf(stderr, "Problems in allocating clovers[%ld] (%s, %d)\n", r, __FILE__, __LINE__);
         exit(EXIT_FAILURE);
         }
       for(i=0; i<STDIM; i++)
          {
          err=posix_memalign((void**)&(GC->clover_array[r][i]), (size_t)DOUBLE_ALIGN, STDIM*sizeof(GAUGE_GROUP));
          if(err!=0)
            {
            fprintf(stderr, "Problems in allocating clovers[%ld][%d] (%s, %d)\n", r, i, __FILE__, __LINE__);
            exit(EXIT_FAILURE);
            }
          }
       }
    }
  }


// free the clovers arrays
void end_clover_array(Gauge_Conf *GC,
                      GParam const * const param)
  {
  int i;
  long r;

  for(r=0; r<param->d_volume; r++)
     {
     for(i=0; i<STDIM; i++)
        {
        free(GC->clover_array[r][i]);
        }
     free(GC->clover_array[r]);
     }
  free(GC->clover_array);
  }


#endif
