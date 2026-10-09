#ifndef GPARAM_C
#define GPARAM_C

#include "../include/macro.h"
#include "../include/endianness.h"
#include "../include/gparam.h"

#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

// skip blanks (spaces, tabs, newlines, the \r of CRLF files) and comments, from # to the end of the
// line, up to the next word or the end of the file
void remove_white_line_and_comments(FILE *input)
{
  int c;

  do
  {
    c = getc(input);
    if (c == '#')
    {
      do
      {
        c = getc(input);
      } while (c != '\n' && c != EOF); // a comment can end the file without a newline
    }
  } while (c != EOF && isspace(c));
  ungetc(c, input); // no-op at EOF
}

// defaults of the optional parameters, and zeros for arrays that the input may fill only in part
static void set_default_parameters(GParam *param)
{
  int i;

  // ml_step[0] = 0 means that the multilevel is not used, and skips its checks
  for (i = 0; i < NLEVELS; i++)
  {
    param->d_ml_step[i] = 0;
  }

  for (i = 0; i < NCOLOR; i++)
  {
    param->d_h[i] = 0.0;
  }
  param->d_theta = 0.0;

  // defect_dir = -1: no defect
  param->d_defect_dir = -1;
  for (i = 0; i < STDIM - 1; i++)
  {
    param->d_L_defect[i] = 0;
  }
  param->d_N_replica_pt = 1;

  // no hierarchical update
  param->d_N_hierarc_levels = 0;
  param->d_L_rect = NULL;
  param->d_N_sweep_rect = NULL;

  param->d_flow_evolutions = 0;
  param->d_flow_between = 0;
  param->d_flow_steps = 0;
  param->d_flow_dmeas = 0;
  param->d_flow_beta_target = 6.0;
  param->d_flow_beta_t_target = 6.0;
  param->d_flow_bc_beta0 = 0.0;
  param->d_flow_protocol_type = 0;
  // allocated by the flow mains only (init_start_end_protocol_*, init_protocol, init_*smearing_parameter)
  param->d_flow_protocol_start = NULL;
  param->d_flow_protocol_end = NULL;
  param->d_flow_protocol = NULL;
  param->d_SNF_rho = NULL;

  // do not measure chi_prime and the time profile of the topological charge
  param->d_chi_prime_meas = 0;
  param->d_topcharge_tprof_meas = 0;
}

// input keywords are matched exactly, so that a misspelled keyword stops the run
static int key_is(char const *key, char const *name)
{
  return strcmp(key, name) == 0;
}

// stop if the value of the keyword key could not be read (err is the return value of fscanf)
static void check_read(int err, char const *key, char const *in_file)
{
  if (err != 1)
  {
    fprintf(stderr, "Error in reading the value of %s in the file %s (%s, %d)\n", key, in_file, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
}

static void read_ints(FILE *input, char const *in_file, char const *key, int *values, int n)
{
  int i;

  for (i = 0; i < n; i++)
  {
    check_read(fscanf(input, "%d", &values[i]), key, in_file);
  }
}

static void read_int(FILE *input, char const *in_file, char const *key, int *value)
{
  read_ints(input, in_file, key, value, 1);
}

// read an integer that has to lie in [min, max]
static void read_int_in_range(FILE *input, char const *in_file, char const *key, int *value, int min, int max)
{
  read_int(input, in_file, key, value);
  if (*value < min || *value > max)
  {
    fprintf(stderr, "Error in reading the file %s: %s must be between %d and %d (%s, %d)\n",
            in_file, key, min, max, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
}

static void read_doubles(FILE *input, char const *in_file, char const *key, double *values, int n)
{
  int i;

  for (i = 0; i < n; i++)
  {
    check_read(fscanf(input, "%lf", &values[i]), key, in_file);
  }
}

static void read_double(FILE *input, char const *in_file, char const *key, double *value)
{
  read_doubles(input, in_file, key, value, 1);
}

#define STRINGIFY_(x) #x
#define STRINGIFY(x) STRINGIFY_(x)

// fscanf(input, "%s", word) into a string of STD_STRING_LENGTH characters (keywords and the file names
// in GParam); stops if the word in the file does not fit. Returns the value of fscanf.
static int scan_word(FILE *input, char const *in_file, char *word)
{
  // one character more than word, to detect the words that do not fit;
  // the format is "%150s" (STD_STRING_LENGTH has to be a plain number)
  char buffer[STD_STRING_LENGTH + 1];
  int const err = fscanf(input, "%" STRINGIFY(STD_STRING_LENGTH) "s", buffer);

  if (err == 1)
  {
    if (strlen(buffer) >= STD_STRING_LENGTH)
    {
      fprintf(stderr, "Error in reading the file %s: %.40s... is longer than %d characters (%s, %d)\n",
              in_file, buffer, STD_STRING_LENGTH - 1, __FILE__, __LINE__);
      exit(EXIT_FAILURE);
    }
    strcpy(word, buffer);
  }
  return err;
}

static void read_string(FILE *input, char const *in_file, char const *key, char *value)
{
  check_read(scan_word(input, in_file, value), key, in_file);
}

// hierarc_upd N  L_rect[0] ... L_rect[N-1]  N_sweep_rect[0] ... N_sweep_rect[N-1]
static void read_hierarc_params(FILE *input, char const *in_file, char const *key, GParam *param)
{
  int n;

  // a repeated hierarc_upd line replaces the previous one
  free(param->d_L_rect);
  free(param->d_N_sweep_rect);
  param->d_L_rect = NULL;
  param->d_N_sweep_rect = NULL;

  read_int(input, in_file, key, &param->d_N_hierarc_levels);
  n = param->d_N_hierarc_levels;
  if (n > 0)
  {
    if (posix_memalign((void **)&(param->d_L_rect), (size_t)INT_ALIGN, (size_t)n * sizeof(int)) != 0
        || posix_memalign((void **)&(param->d_N_sweep_rect), (size_t)INT_ALIGN, (size_t)n * sizeof(int)) != 0)
    {
      fprintf(stderr, "Problems in allocating hierarchical update parameters! (%s, %d)\n", __FILE__, __LINE__);
      exit(EXIT_FAILURE);
    }
    read_ints(input, in_file, key, param->d_L_rect, n);
    read_ints(input, in_file, key, param->d_N_sweep_rect, n);
  }
}

// read the value(s) that follow the keyword key in the input file
static void read_keyword_value(FILE *input, char const *in_file, char const *key, GParam *param)
{
  // lattice and action
  if (key_is(key, "size"))
    read_ints(input, in_file, key, param->d_size, STDIM);
  else if (key_is(key, "beta"))
    read_double(input, in_file, key, &param->d_beta);
  else if (key_is(key, "beta_t"))
    read_double(input, in_file, key, &param->d_beta_t);
  else if (key_is(key, "anisotropic"))
    read_int(input, in_file, key, &param->d_anisotropic);
  else if (key_is(key, "htracedef"))
    read_doubles(input, in_file, key, param->d_h, NCOLOR / 2);
  else if (key_is(key, "theta"))
    read_double(input, in_file, key, &param->d_theta);

  // Monte Carlo
  else if (key_is(key, "sample"))
    read_int(input, in_file, key, &param->d_sample);
  else if (key_is(key, "thermal"))
    read_int(input, in_file, key, &param->d_thermal);
  else if (key_is(key, "overrelax"))
    read_int(input, in_file, key, &param->d_overrelax);
  else if (key_is(key, "measevery"))
    read_int(input, in_file, key, &param->d_measevery);
  else if (key_is(key, "start"))
    read_int(input, in_file, key, &param->d_start);
  else if (key_is(key, "saveconf_back_every"))
    read_int(input, in_file, key, &param->d_saveconf_back_every);
  else if (key_is(key, "saveconf_analysis_every"))
    read_int(input, in_file, key, &param->d_saveconf_analysis_every);
  else if (key_is(key, "epsilon_metro"))
    read_double(input, in_file, key, &param->d_epsilon_metro);
  else if (key_is(key, "randseed"))
    check_read(fscanf(input, "%u", &param->d_randseed), key, in_file);

  // measurements: cooling, topological charge, gradient flow
  else if (key_is(key, "coolsteps"))
    read_int(input, in_file, key, &param->d_coolsteps);
  else if (key_is(key, "coolrepeat"))
    read_int(input, in_file, key, &param->d_coolrepeat);
  else if (key_is(key, "chi_prime_meas"))
    read_int_in_range(input, in_file, key, &param->d_chi_prime_meas, 0, 1);
  else if (key_is(key, "topcharge_tprof_meas"))
    read_int_in_range(input, in_file, key, &param->d_topcharge_tprof_meas, 0, 1);
  else if (key_is(key, "gfstep")) // integration step
    read_double(input, in_file, key, &param->d_gfstep);
  else if (key_is(key, "num_gfsteps")) // number of integration steps
    read_int(input, in_file, key, &param->d_ngfsteps);
  else if (key_is(key, "gf_meas_each"))
    read_int(input, in_file, key, &param->d_gf_meas_each);

  // multilevel
  else if (key_is(key, "multihit"))
    read_int(input, in_file, key, &param->d_multihit);
  else if (key_is(key, "ml_step"))
    read_ints(input, in_file, key, param->d_ml_step, NLEVELS);
  else if (key_is(key, "ml_upd"))
    read_ints(input, in_file, key, param->d_ml_upd, NLEVELS);
  else if (key_is(key, "ml_level0_repeat"))
    read_int(input, in_file, key, &param->d_ml_level0_repeat);
  else if (key_is(key, "dist_poly"))
    read_int(input, in_file, key, &param->d_dist_poly);
  else if (key_is(key, "transv_dist"))
    read_int(input, in_file, key, &param->d_trasv_dist);
  else if (key_is(key, "plaq_dir"))
    read_ints(input, in_file, key, param->d_plaq_dir, 2);

  // defect (PTBC) and hierarchical update
  else if (key_is(key, "defect_dir"))
    read_int_in_range(input, in_file, key, &param->d_defect_dir, 0, STDIM - 1);
  else if (key_is(key, "defect_size"))
    read_ints(input, in_file, key, param->d_L_defect, STDIM - 1);
  else if (key_is(key, "hierarc_upd"))
    read_hierarc_params(input, in_file, key, param);

  // non-equilibrium evolutions (Jarzynski, SNF)
  else if (key_is(key, "flow_beta_target"))
    read_double(input, in_file, key, &param->d_flow_beta_target);
  else if (key_is(key, "flow_beta_t_target"))
    read_double(input, in_file, key, &param->d_flow_beta_t_target);
  else if (key_is(key, "flow_bc_beta0"))
    read_double(input, in_file, key, &param->d_flow_bc_beta0);
  else if (key_is(key, "num_flow_ev"))
    read_int(input, in_file, key, &param->d_flow_evolutions);
  else if (key_is(key, "num_flow_steps"))
    read_int(input, in_file, key, &param->d_flow_steps);
  else if (key_is(key, "num_flow_between"))
    read_int(input, in_file, key, &param->d_flow_between);
  else if (key_is(key, "num_flow_dmeas"))
    read_int(input, in_file, key, &param->d_flow_dmeas);
  else if (key_is(key, "protocol_type"))
    read_int(input, in_file, key, &param->d_flow_protocol_type);

  // multicanonic
  else if (key_is(key, "grid_step"))
    read_double(input, in_file, key, &param->d_grid_step);
  else if (key_is(key, "grid_max"))
    read_double(input, in_file, key, &param->d_grid_max);

  // file names
  else if (key_is(key, "conf_file"))
    read_string(input, in_file, key, param->d_conf_file);
  else if (key_is(key, "data_file"))
    read_string(input, in_file, key, param->d_data_file);
  else if (key_is(key, "work_file"))
    read_string(input, in_file, key, param->d_work_file);
  else if (key_is(key, "log_file"))
    read_string(input, in_file, key, param->d_log_file);
  else if (key_is(key, "protocol_file"))
    read_string(input, in_file, key, param->d_protocol_file);
  else if (key_is(key, "smearingrho_file"))
    read_string(input, in_file, key, param->d_smearingrho_file);
  else if (key_is(key, "chiprime_data_file"))
    read_string(input, in_file, key, param->d_chiprime_file);
  else if (key_is(key, "topcharge_tprof_file"))
    read_string(input, in_file, key, param->d_topcharge_tprof_file);
  else if (key_is(key, "ml_file"))
    read_string(input, in_file, key, param->d_ml_file);
  else if (key_is(key, "swap_acc_file"))
    read_string(input, in_file, key, param->d_swap_acc_file);
  else if (key_is(key, "swap_track_file"))
    read_string(input, in_file, key, param->d_swap_tracking_file);
  else if (key_is(key, "multicanonic_acc_file"))
    read_string(input, in_file, key, param->d_multicanonic_acc_file);
  else if (key_is(key, "topo_potential_file"))
    read_string(input, in_file, key, param->d_topo_potential_file);

  else
  {
    fprintf(stderr, "Error: unrecognized option %s in the file %s (%s, %d)\n", key, in_file, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
}

static int end_of_file(FILE *input)
{
  int const c = getc(input);

  if (c == EOF)
  {
    return 1;
  }
  ungetc(c, input);
  return 0;
}

// multilevel: ml_step[0] divides size[0], each ml_step[i] divides ml_step[i-1] and is smaller, and all are > 1
static void check_multilevel_steps(GParam const *param)
{
  int i;

  if (param->d_ml_step[0] == 0) // multilevel not used
  {
    return;
  }

  if (param->d_size[0] % param->d_ml_step[0] || param->d_size[0] < param->d_ml_step[0])
  {
    fprintf(stderr, "Error: size[0] has to be divisible by ml_step[0] and satisfy ml_step[0]<=size[0] (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  for (i = 1; i < NLEVELS; i++)
  {
    if (param->d_ml_step[i - 1] % param->d_ml_step[i] || param->d_ml_step[i - 1] <= param->d_ml_step[i])
    {
      fprintf(stderr, "Error: ml_step[%d] has to be divisible by ml_step[%d] and larger than it (%s, %d)\n", i - 1, i, __FILE__, __LINE__);
      exit(EXIT_FAILURE);
    }
  }
  if (param->d_ml_step[NLEVELS - 1] == 1)
  {
    fprintf(stderr, "Error: ml_step[%d] has to be larger than 1 (%s, %d)\n", NLEVELS - 1, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
}

static void check_lattice_sizes(GParam const *param)
{
  int i;

#ifdef OPENMP_MODE
  // the even/odd parallel updates need even sides
  for (i = 0; i < STDIM; i++)
  {
    if (param->d_size[i] % 2 != 0)
    {
      fprintf(stderr, "Error: size[%d] is not even.\n", i);
      fprintf(stderr, "When using OpenMP all the sides of the lattice have to be even! (%s, %d)\n", __FILE__, __LINE__);
      exit(EXIT_FAILURE);
    }
  }
#endif

  for (i = 0; i < STDIM; i++)
  {
    if (param->d_size[i] == 1)
    {
      fprintf(stderr, "Error: all sizes has to be larger than 1: the totally reduced case is not implemented! (%s, %d)\n", __FILE__, __LINE__);
      exit(EXIT_FAILURE);
    }
  }
}

// L_defect[k] is the size of the defect along the k-th direction orthogonal to defect_dir, in increasing
// order (as perp_dir in geometry.c): for defect_dir = 1, along t, y, z
static void check_defect_size(GParam const *param)
{
  int dir, k;

  if (param->d_defect_dir < 0) // no defect
  {
    for (k = 0; k < STDIM - 1; k++)
    {
      if (param->d_L_defect[k] != 0)
      {
        fprintf(stderr, "Error: defect_size needs defect_dir (%s, %d)\n", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
      }
    }
    return;
  }

  k = 0;
  for (dir = 0; dir < STDIM; dir++)
  {
    if (dir == param->d_defect_dir)
    {
      continue;
    }
    if (param->d_L_defect[k] > param->d_size[dir])
    {
      fprintf(stderr, "Error: defect_size[%d] = %d is larger than size[%d] = %d (defect_dir %d) (%s, %d)\n",
              k, param->d_L_defect[k], dir, param->d_size[dir], param->d_defect_dir, __FILE__, __LINE__);
      exit(EXIT_FAILURE);
    }
    k++;
  }
}

// the input file is a sequence of keywords, each followed by its value(s);
// empty lines and comments (from # to the end of the line) are skipped
void readinput(char *in_file, GParam *param)
{
  FILE *input;
  char key[STD_STRING_LENGTH];
  int err;

  set_default_parameters(param);

  input = fopen(in_file, "r");
  if (input == NULL)
  {
    fprintf(stderr, "Error in opening the file %s (%s, %d)\n", in_file, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }

  do
  {
    remove_white_line_and_comments(input);

    err = scan_word(input, in_file, key);
    if (err != 1)
    {
      fprintf(stderr, "Error in reading the file %s, err=%d (%s, %d)\n", in_file, err, __FILE__, __LINE__);
      exit(EXIT_FAILURE);
    }
    read_keyword_value(input, in_file, key, param);

    remove_white_line_and_comments(input);
  } while (!end_of_file(input));

  fclose(input);

  check_multilevel_steps(param);
  check_lattice_sizes(param);
  check_defect_size(param);

  init_derived_constants(param);
}

void init_derived_constants(GParam *param)
{
  int i;

  // derived constants
  param->d_volume = 1;
  for (i = 0; i < STDIM; i++)
  {
    (param->d_volume) *= (param->d_size[i]);
  }

  param->d_space_vol = 1;
  // direction 0 is time
  for (i = 1; i < STDIM; i++)
  {
    (param->d_space_vol) *= (param->d_size[i]);
  }

  param->d_inv_vol = 1.0 / ((double)param->d_volume);
  param->d_inv_space_vol = 1.0 / ((double)param->d_space_vol);

  // volume of the defect
  param->d_volume_defect = 1;
  for (i = 0; i < STDIM - 1; i++)
  {
    param->d_volume_defect *= param->d_L_defect[i];
  }

  // number of grid points (multicanonic only)
  param->d_n_grid = (int)((2.0 * param->d_grid_max / param->d_grid_step) + 1.0);

  
  // for isotropic simulations only beta_t is set to be beta
  if (param->d_anisotropic == 0)
  {
    param->d_beta_t = param->d_beta;
    param->d_flow_beta_t_target = param->d_flow_beta_target;
  }
}

static void check_flow_steps(GParam const *param)
{
  if (param->d_flow_steps < 1)
  {
    fprintf(stderr, "Error: num_flow_steps = %d, it has to be at least 1 (%s, %d)\n",
            param->d_flow_steps, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
}

// checks of the input of the flows in beta, to be called after readinput
void check_flow_beta_input(GParam const *param)
{
  check_flow_steps(param);

  // the flows in beta measure every num_flow_dmeas steps of an evolution
  if (param->d_flow_dmeas < 1)
  {
    fprintf(stderr, "Error: num_flow_dmeas = %d, it has to be at least 1 (%s, %d)\n",
            param->d_flow_dmeas, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
}

// checks of the input of the flows in the boundary conditions, to be called after readinput
void check_flow_bc_input(GParam const *param)
{
  int k;

  check_flow_steps(param);

  // the flows in the boundary conditions need a defect (defect_dir = -1, the default, means no defect)
  if (param->d_defect_dir < 0)
  {
    fprintf(stderr, "Error: a flow in the boundary conditions needs defect_dir and defect_size (%s, %d)\n",
            __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  // and a non-empty one (defect_size defaults to 0)
  for (k = 0; k < STDIM - 1; k++)
  {
    if (param->d_L_defect[k] < 1)
    {
      fprintf(stderr, "Error: defect_size[%d] = %d, it has to be at least 1 (%s, %d)\n",
              k, param->d_L_defect[k], __FILE__, __LINE__);
      exit(EXIT_FAILURE);
    }
  }
}

// allocate the starting and ending values of the npar protocol parameters
static void alloc_start_end_protocol(GParam *param, int npar)
{
  int err;

  err = posix_memalign((void **)&(param->d_flow_protocol_start), (size_t)DOUBLE_ALIGN, (size_t)npar * sizeof(double));
  if (err != 0)
  {
    fprintf(stderr, "Problems in allocating protocol parameters! (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  err = posix_memalign((void **)&(param->d_flow_protocol_end), (size_t)DOUBLE_ALIGN, (size_t)npar * sizeof(double));
  if (err != 0)
  {
    fprintf(stderr, "Problems in allocating protocol parameters! (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
}

void init_start_end_protocol_beta(GParam *param, int npar)
{
  alloc_start_end_protocol(param, npar);

  if (param->d_anisotropic)
  {
    param->d_flow_protocol_start[0] = param->d_beta;
    param->d_flow_protocol_start[1] = param->d_beta_t;
    param->d_flow_protocol_end[0] = param->d_flow_beta_target;
    param->d_flow_protocol_end[1] = param->d_flow_beta_t_target;
  }
  else
  {
    param->d_flow_protocol_start[0] = param->d_beta;
    param->d_flow_protocol_end[0] = param->d_flow_beta_target;
  }
}

void init_start_end_protocol_bc(GParam *param)
{
  alloc_start_end_protocol(param, 1);

  param->d_flow_protocol_start[0] = param->d_flow_bc_beta0;
  param->d_flow_protocol_end[0] = 1.0;
}

// d_flow_protocol[p * d_flow_steps + i] = value of the parameter p after step i (p = 0: beta, p = 1: beta_t)
void init_protocol(GParam *param, int npar)
{
  FILE *input_protocol;
  double temp_d;
  int i, p;
  int err;

  err = posix_memalign((void **)&(param->d_flow_protocol), (size_t)DOUBLE_ALIGN, (size_t)param->d_flow_steps * npar * sizeof(double));
  if (err != 0)
  {
    fprintf(stderr, "Problems in allocating protocol parameters! (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }

  if (param->d_flow_protocol_type)
  {
    input_protocol = fopen(param->d_protocol_file, "r"); // open the input protocol file

    if (input_protocol == NULL)
    {
      fprintf(stderr, "Error in opening the file %s (%s, %d)\n", param->d_protocol_file, __FILE__, __LINE__);
      exit(EXIT_FAILURE);
    }
    else
    {
      for (p = 0; p < npar; p++)
        for (i = 0; i < param->d_flow_steps; i++)
        {
          err = fscanf(input_protocol, "%lf", &temp_d);
          if (err != 1)
          {
            fprintf(stderr, "Error in reading the file %s (%s, %d)\n", param->d_protocol_file, __FILE__, __LINE__);
            exit(EXIT_FAILURE);
          }
          param->d_flow_protocol[p * param->d_flow_steps + i] = temp_d;
        }
      fclose(input_protocol);
    }
  }
  else
  {
    for (p = 0; p < npar; p++)
      for (i = 0; i < param->d_flow_steps; i++)
        param->d_flow_protocol[p * param->d_flow_steps + i] = (double)((param->d_flow_protocol_end[p] - param->d_flow_protocol_start[p]) * ((double)(i + 1)) / param->d_flow_steps + param->d_flow_protocol_start[p]);
  }
}

// flows in beta: set beta (and beta_t, if anisotropic) to their values after step `step` of the protocol
void set_flow_beta(GParam *param, int step)
{
  param->d_beta = param->d_flow_protocol[step];
  if (param->d_anisotropic != 0)
  {
    param->d_beta_t = param->d_flow_protocol[param->d_flow_steps + step];
  }
}

void init_smearing_parameter(GParam *param)
{
  FILE *input_smearingrho;
  double temp_d;
  int i;
  int err;

  err = posix_memalign((void **)&(param->d_SNF_rho), (size_t)DOUBLE_ALIGN, (size_t)param->d_flow_steps * sizeof(double));
  if (err != 0)
  {
    fprintf(stderr, "Problems in allocating protocol parameters! (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }

  input_smearingrho = fopen(param->d_smearingrho_file, "r"); // open the input smearing rho file

  if (input_smearingrho == NULL)
  {
    fprintf(stderr, "Error in opening the file %s (%s, %d)\n", param->d_smearingrho_file, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  else
  {
    for (i = 0; i < param->d_flow_steps; i++)
    {
      err = fscanf(input_smearingrho, "%lf", &temp_d);
      if (err != 1)
      {
        fprintf(stderr, "Error in reading the file %s (%s, %d)\n", param->d_smearingrho_file, __FILE__, __LINE__);
        exit(EXIT_FAILURE);
      }
      param->d_SNF_rho[i] = temp_d;
    }
    fclose(input_smearingrho);
  }
}

// d_SNF_rho[((i * STDIM + mu) * rect_vol + s) * 2(STDIM-1) + p] = rho of step i, link direction mu, site s of
// the smearing rectangle (rect_sites order), staple p (calcstaples_wilson_nosum order)
void init_defect_smearing_parameter(GParam *param, long rect_vol)
{
  FILE *input_smearingrho;
  double temp_d;
  int i, mu, s, p;
  long err;

  err = posix_memalign((void **)&(param->d_SNF_rho), (size_t)DOUBLE_ALIGN, (size_t)param->d_flow_steps * 2 * (STDIM - 1) * STDIM * rect_vol * sizeof(double));
  if (err != 0)
  {
    fprintf(stderr, "Problems in allocating protocol parameters! (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }

  input_smearingrho = fopen(param->d_smearingrho_file, "rb"); // open the input smearing rho binary file

  if (input_smearingrho == NULL)
  {
    fprintf(stderr, "Error in opening the file %s (%s, %d)\n", param->d_smearingrho_file, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  else
  {
    // int endianness = endian();
    for (i = 0; i < param->d_flow_steps; i++)
      for (mu = 0; mu < STDIM; mu++)
        for (s = 0; s < rect_vol; s++)
          for (p = 0; p < 2 * (STDIM - 1); p++)
          {
            long rho_index = 2 * (STDIM - 1) * rect_vol * STDIM * i + 2 * (STDIM - 1) * rect_vol * mu + 2 * (STDIM - 1) * s + p;
            err = fread(&temp_d, sizeof(double), 1, input_smearingrho);
            if (err != 1)
            {
              fprintf(stderr, "Error in reading the file %s (%s, %d)\n", param->d_smearingrho_file, __FILE__, __LINE__);
              exit(EXIT_FAILURE);
            }

            // if (endianness == 0)
            //   SwapBytesDouble(&temp_d);

            param->d_SNF_rho[rho_index] = temp_d;
          }
    fclose(input_smearingrho);
  }
}

// rho of the defect smearing at step `step`: STDIM * rect_vol * 2(STDIM-1) values, ordered as in
// init_defect_smearing_parameter (the layout defect_stout_smearing_update expects)
double *defect_smearing_rho(GParam const *param, long rect_vol, int step)
{
  return param->d_SNF_rho + 2 * (STDIM - 1) * rect_vol * STDIM * step;
}

// initialize data file
void init_data_file(FILE **dataf, FILE **chiprimef, FILE **topchar_tprof_f, GParam const *const param)
{
  int i;

  if (param->d_start == 2)
  {
    // open std data file (plaquette, polyakov, topological charge)
    *dataf = fopen(param->d_data_file, "r");
    if (*dataf != NULL) // file exists
    {
      fclose(*dataf);
      *dataf = fopen(param->d_data_file, "a");
    }
    else
    {
      *dataf = fopen(param->d_data_file, "w");
      fprintf(*dataf, "%d ", STDIM);
      for (i = 0; i < STDIM; i++)
        fprintf(*dataf, "%d ", param->d_size[i]);
      fprintf(*dataf, "\n");
    }
    // open chi prime data file
    if (param->d_chi_prime_meas == 1)
    {
      *chiprimef = fopen(param->d_chiprime_file, "r");
      if (*chiprimef != NULL) // file exists
      {
        fclose(*chiprimef);
        *chiprimef = fopen(param->d_chiprime_file, "a");
      }
      else
      {
        *chiprimef = fopen(param->d_chiprime_file, "w");
        fprintf(*chiprimef, "%d ", STDIM);
        for (i = 0; i < STDIM; i++)
          fprintf(*chiprimef, "%d ", param->d_size[i]);
        fprintf(*chiprimef, "\n");
      }
    }
    else
    {
      (void)chiprimef;
    }
    // open topocharge_tprof data file
    if (param->d_topcharge_tprof_meas == 1)
    {
      *topchar_tprof_f = fopen(param->d_topcharge_tprof_file, "r");
      if (*topchar_tprof_f != NULL) // file exists
      {
        fclose(*topchar_tprof_f);
        *topchar_tprof_f = fopen(param->d_topcharge_tprof_file, "a");
      }
      else
      {
        *topchar_tprof_f = fopen(param->d_topcharge_tprof_file, "w");
        fprintf(*topchar_tprof_f, "%d ", STDIM);
        for (i = 0; i < STDIM; i++)
          fprintf(*topchar_tprof_f, "%d ", param->d_size[i]);
        fprintf(*topchar_tprof_f, "\n");
      }
    }
    else
    {
      (void)topchar_tprof_f;
    }
  }
  else
  {
    // open std data file
    *dataf = fopen(param->d_data_file, "w");
    fprintf(*dataf, "%d ", STDIM);
    for (i = 0; i < STDIM; i++)
    {
      fprintf(*dataf, "%d ", param->d_size[i]);
    }
    fprintf(*dataf, "\n");
    // open chi prime data file
    if (param->d_chi_prime_meas == 1)
    {
      *chiprimef = fopen(param->d_chiprime_file, "w");
      fprintf(*chiprimef, "%d ", STDIM);
      for (i = 0; i < STDIM; i++)
        fprintf(*chiprimef, "%d ", param->d_size[i]);
      fprintf(*chiprimef, "\n");
    }
    else
    {
      (void)chiprimef;
    }
    // open topocharge_tprof data file
    if (param->d_topcharge_tprof_meas == 1)
    {
      *topchar_tprof_f = fopen(param->d_topcharge_tprof_file, "w");
      fprintf(*topchar_tprof_f, "%d ", STDIM);
      for (i = 0; i < STDIM; i++)
        fprintf(*topchar_tprof_f, "%d ", param->d_size[i]);
      fprintf(*topchar_tprof_f, "\n");
    }
    else
    {
      (void)topchar_tprof_f;
    }
  }
  fflush(*dataf);
  if (param->d_chi_prime_meas == 1)
    fflush(*chiprimef);
  else
  {
    (void)chiprimef;
  }
  if (param->d_topcharge_tprof_meas == 1)
    fflush(*topchar_tprof_f);
  else
  {
    (void)topchar_tprof_f;
  }
}

void init_work_file(FILE **workfilep, GParam const *const param)
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
    if (param->d_anisotropic != 0)
      fprintf(*workfilep, " %f ", param->d_beta_t);
    fprintf(*workfilep, "%d ", STDIM);
    for (i = 0; i < STDIM; i++)
      fprintf(*workfilep, "%d ", param->d_size[i]);
    fprintf(*workfilep, "\n");
  }
  fflush(*workfilep);
}

// free allocated memory for hierarc update parameters
void free_hierarc_params(GParam *param)
{
  if (param->d_N_hierarc_levels == 0)
  {
    (void)param; // to avoid compiler warning about unused variable
  }
  else
  {
    free(param->d_L_rect);
    free(param->d_N_sweep_rect);
  }
}

// free the protocol and smearing parameters of the flows (NULL if not allocated)
void free_flow_params(GParam *param)
{
  free(param->d_flow_protocol_start);
  free(param->d_flow_protocol_end);
  free(param->d_flow_protocol);
  free(param->d_SNF_rho);
  param->d_flow_protocol_start = NULL;
  param->d_flow_protocol_end = NULL;
  param->d_flow_protocol = NULL;
  param->d_SNF_rho = NULL;
}

// print simulation parameters

void print_parameters_local(GParam const *const param, time_t time_start, time_t time_end)
{
  FILE *fp;
  int i;
  double diff_sec;

  fp = fopen(param->d_log_file, "w");
  fprintf(fp, "+-----------------------------------------+\n");
  fprintf(fp, "| Simulation details for yang_mills_local |\n");
  fprintf(fp, "+-----------------------------------------+\n\n");

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

  fprintf(fp, "anisotropy: %d\n", param->d_anisotropic);
  fprintf(fp, "beta: %.10lf\n", param->d_beta);
  if (param->d_anisotropic != 0)
    fprintf(fp, "beta_t: %.10lf\n", param->d_beta_t);
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

// void print_parameters_local_pt_multicanonic(GParam const * const param, time_t time_start, time_t time_end)
//     {
//     FILE *fp;
//     int i;
//     double diff_sec;

//     fp=fopen(param->d_log_file, "w");
//     fprintf(fp, "+---------------------------------------------------------+\n");
//     fprintf(fp, "| Simulation details for yang_mills_local_pt_multicanonic |\n");
//     fprintf(fp, "+---------------------------------------------------------+\n\n");

//     #ifdef OPENMP_MODE
//      fprintf(fp, "using OpenMP with %d threads\n\n", NTHREADS);
//     #endif

//     fprintf(fp, "number of colors: %d\n", NCOLOR);
//     fprintf(fp, "spacetime dimensionality: %d\n\n", STDIM);

//     fprintf(fp, "lattice: %d", param->d_size[0]);
//     for(i=1; i<STDIM; i++)
//        {
//        fprintf(fp, "x%d", param->d_size[i]);
//        }
//     fprintf(fp, "\n\n");

// 	fprintf(fp, "defect dir: %d\n", param->d_defect_dir);
//     fprintf(fp, "defect: %d", param->d_L_defect[0]);
//     for(i=1; i<STDIM-1; i++)
//        {
//        fprintf(fp, "x%d", param->d_L_defect[i]);
//        }
//     fprintf(fp, "\n\n");
//     fprintf(fp,"number of copies used in parallel tempering: %d\n", param->d_N_replica_pt);
// 		fprintf(fp,"boundary condition constants: ");
// 		for(i=0;i<param->d_N_replica_pt;i++)
// 			fprintf(fp,"%lf ",param->d_pt_bound_cond_coeff[i]);
// 		fprintf(fp,"\n");
// 		fprintf(fp,"number of hierarchical levels: %d\n", param->d_N_hierarc_levels);
// 		if(param->d_N_hierarc_levels>0)
// 			{
// 			fprintf(fp,"extention of rectangles: ");
// 			for(i=0;i<param->d_N_hierarc_levels;i++)
// 				{
// 				fprintf(fp,"%d ", param->d_L_rect[i]);
// 				}
// 			fprintf(fp,"\n");
// 			fprintf(fp,"number of sweeps per hierarchical level: ");
// 			for(i=0;i<param->d_N_hierarc_levels;i++)
// 				{
// 				fprintf(fp,"%d ", param->d_N_sweep_rect[i]);
// 				}
// 			}
// 		fprintf(fp,"\n\n");

// 		fprintf(fp,"Multicanonic topo-potential read from file %s\nPotential defined on a grid with step=%.10lf and max=%.10lf\n", param->d_topo_potential_file, param->d_grid_step, param->d_grid_max);

// 		fprintf(fp,"\n\n");

//     fprintf(fp, "beta: %.10lf\n", param->d_beta);
//     fprintf(fp, "beta_t: %.10lf\n", param->d_beta_t);
//     #ifdef THETA_MODE
//       fprintf(fp, "theta: %.10lf\n", param->d_theta);
//     #endif
//     fprintf(fp, "\n");

//     fprintf(fp, "sample:    %d\n", param->d_sample);
//     fprintf(fp, "thermal:   %d\n", param->d_thermal);
//     fprintf(fp, "overrelax: %d\n", param->d_overrelax);
//     fprintf(fp, "measevery: %d\n", param->d_measevery);
//     fprintf(fp, "\n");

//     fprintf(fp, "start:                   %d\n", param->d_start);
//     fprintf(fp, "saveconf_back_every:     %d\n", param->d_saveconf_back_every);
//     fprintf(fp, "saveconf_analysis_every: %d\n", param->d_saveconf_analysis_every);
//     fprintf(fp, "\n");

//     fprintf(fp, "coolsteps:      %d\n", param->d_coolsteps);
//     fprintf(fp, "coolrepeat:     %d\n", param->d_coolrepeat);
//     fprintf(fp, "\n");

//     fprintf(fp, "randseed: %u\n", param->d_randseed);
//     fprintf(fp, "\n");

//     diff_sec = difftime(time_end, time_start);
//     fprintf(fp, "Simulation time: %.3lf seconds\n", diff_sec );
//     fprintf(fp, "\n");

//     if(endian()==0)
//       {
//       fprintf(fp, "Little endian machine\n\n");
//       }
//     else
//       {
//       fprintf(fp, "Big endian machine\n\n");
//       }

//     fclose(fp);
//     }

// void print_parameters_local_pt(GParam const * const param, time_t time_start, time_t time_end)
//     {
//     FILE *fp;
//     int i;
//     double diff_sec;

//     fp=fopen(param->d_log_file, "w");
//     fprintf(fp, "+--------------------------------------------+\n");
//     fprintf(fp, "| Simulation details for yang_mills_local_pt |\n");
//     fprintf(fp, "+--------------------------------------------+\n\n");

//     #ifdef OPENMP_MODE
//      fprintf(fp, "using OpenMP with %d threads\n\n", NTHREADS);
//     #endif

//     fprintf(fp, "number of colors: %d\n", NCOLOR);
//     fprintf(fp, "spacetime dimensionality: %d\n\n", STDIM);

//     fprintf(fp, "lattice: %d", param->d_size[0]);
//     for(i=1; i<STDIM; i++)
//        {
//        fprintf(fp, "x%d", param->d_size[i]);
//        }
//     fprintf(fp, "\n\n");

// 	fprintf(fp, "defect dir: %d\n", param->d_defect_dir);
//     fprintf(fp, "defect: %d", param->d_L_defect[0]);
//     for(i=1; i<STDIM-1; i++)
//        {
//        fprintf(fp, "x%d", param->d_L_defect[i]);
//        }
//     fprintf(fp, "\n\n");
//     fprintf(fp,"number of copies used in parallel tempering: %d\n", param->d_N_replica_pt);
// 		fprintf(fp,"boundary condition constants: ");
// 		for(i=0;i<param->d_N_replica_pt;i++)
// 			fprintf(fp,"%lf ",param->d_pt_bound_cond_coeff[i]);
// 		fprintf(fp,"\n");
// 		fprintf(fp,"number of hierarchical levels: %d\n", param->d_N_hierarc_levels);
// 		if(param->d_N_hierarc_levels>0)
// 			{
// 			fprintf(fp,"extention of rectangles: ");
// 			for(i=0;i<param->d_N_hierarc_levels;i++)
// 				{
// 				fprintf(fp,"%d ", param->d_L_rect[i]);
// 				}
// 			fprintf(fp,"\n");
// 			fprintf(fp,"number of sweeps per hierarchical level: ");
// 			for(i=0;i<param->d_N_hierarc_levels;i++)
// 				{
// 				fprintf(fp,"%d ", param->d_N_sweep_rect[i]);
// 				}
// 			}
// 		fprintf(fp,"\n\n");

//     fprintf(fp, "beta: %.10lf\n", param->d_beta);
//     fprintf(fp, "beta_t: %.10lf\n", param->d_beta_t);
//     #ifdef THETA_MODE
//       fprintf(fp, "theta: %.10lf\n", param->d_theta);
//     #endif
//     fprintf(fp, "\n");

//     fprintf(fp, "sample:    %d\n", param->d_sample);
//     fprintf(fp, "thermal:   %d\n", param->d_thermal);
//     fprintf(fp, "overrelax: %d\n", param->d_overrelax);
//     fprintf(fp, "measevery: %d\n", param->d_measevery);
//     fprintf(fp, "\n");

//     fprintf(fp, "start:                   %d\n", param->d_start);
//     fprintf(fp, "saveconf_back_every:     %d\n", param->d_saveconf_back_every);
//     fprintf(fp, "saveconf_analysis_every: %d\n", param->d_saveconf_analysis_every);
//     fprintf(fp, "\n");

//     fprintf(fp, "coolsteps:      %d\n", param->d_coolsteps);
//     fprintf(fp, "coolrepeat:     %d\n", param->d_coolrepeat);
//     fprintf(fp, "\n");

//     fprintf(fp, "randseed: %u\n", param->d_randseed);
//     fprintf(fp, "\n");

//     diff_sec = difftime(time_end, time_start);
//     fprintf(fp, "Simulation time: %.3lf seconds\n", diff_sec );
//     fprintf(fp, "\n");

//     if(endian()==0)
//       {
//       fprintf(fp, "Little endian machine\n\n");
//       }
//     else
//       {
//       fprintf(fp, "Big endian machine\n\n");
//       }

//     fclose(fp);
//     }

void print_parameters_local_flow_bc(GParam const *const param, time_t time_start, time_t time_end)
{
  FILE *fp;
  int i;
  double diff_sec;

  fp = fopen(param->d_log_file, "w");
  fprintf(fp, "+---------------------------------------------------+\n");
  fprintf(fp, "| Simulation details for yang_mills_local_jarzynski/snf_bc |\n");
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
  fprintf(fp, "number of out-of-equilibrium evolutions: %d\n", param->d_flow_evolutions);
  fprintf(fp, "number of out-of-equilibrium steps in each evolution: %d\n", param->d_flow_steps);
  fprintf(fp, "number of relax updates between evolutions: %d\n", param->d_flow_between);
  fprintf(fp, "number of steps between measurements during evolution: %d\n", param->d_flow_dmeas);
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

  fprintf(fp, "anisotropy: %d\n", param->d_anisotropic);
  fprintf(fp, "beta: %.10lf\n", param->d_beta);
  if (param->d_anisotropic != 0)
    fprintf(fp, "beta_t: %.10lf\n", param->d_beta_t);
#ifdef THETA_MODE
  fprintf(fp, "theta: %.10lf\n", param->d_theta);
#endif
  fprintf(fp, "\n");

  fprintf(fp, "thermal:   %d\n", param->d_thermal);
  fprintf(fp, "overrelax: %d\n", param->d_overrelax);
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

void print_parameters_local_flow_beta(GParam const *const param, time_t time_start, time_t time_end)
{
  FILE *fp;
  int i;
  double diff_sec;

  fp = fopen(param->d_log_file, "w");
  fprintf(fp, "+---------------------------------------------------+\n");
  fprintf(fp, "| Simulation details for yang_mills_local_jarzynski/snf_beta |\n");
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

  fprintf(fp, "number of out-of-equilibrium evolutions: %d\n", param->d_flow_evolutions);
  fprintf(fp, "number of out-of-equilibrium steps in each evolution: %d\n", param->d_flow_steps);
  fprintf(fp, "number of relax updates between evolutions: %d\n", param->d_flow_between);
  fprintf(fp, "number of steps between measurements during evolution: %d\n", param->d_flow_dmeas);
  fprintf(fp, "\n\n");

  fprintf(fp, "anisotropy: %d\n", param->d_anisotropic);
  fprintf(fp, "beta_0: %.10lf\n", param->d_beta);
  fprintf(fp, "beta_target: %.10lf\n", param->d_flow_beta_target);
  if (param->d_anisotropic != 0)
    fprintf(fp, "beta_t_0: %.10lf\n", param->d_beta_t);
    fprintf(fp, "beta_t_target: %.10lf\n", param->d_flow_beta_t_target);
#ifdef THETA_MODE
  fprintf(fp, "theta: %.10lf\n", param->d_theta);
#endif
  fprintf(fp, "\n");

  fprintf(fp, "thermal:   %d\n", param->d_thermal);
  fprintf(fp, "overrelax: %d\n", param->d_overrelax);
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
void print_parameters_polycorr_long(GParam *param, time_t time_start, time_t time_end)
{
  FILE *fp;
  int i;
  double diff_sec;

  fp = fopen(param->d_log_file, "w");
  fprintf(fp, "+-------------------------------------------------+\n");
  fprintf(fp, "| Simulation details for yang_mills_polycorr_long |\n");
  fprintf(fp, "+-------------------------------------------------+\n\n");

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

  fprintf(fp, "anisotropy: %d\n", param->d_anisotropic);
  fprintf(fp, "beta: %.10lf\n", param->d_beta);
  if (param->d_anisotropic != 0)
    fprintf(fp, "beta_t: %.10lf\n", param->d_beta_t);
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
  for (i = 0; i < NLEVELS; i++)
  {
    fprintf(fp, "%d ", param->d_ml_step[i]);
  }
  fprintf(fp, "\n");
  fprintf(fp, "updates for levels: ");
  for (i = 0; i < NLEVELS; i++)
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
void print_parameters_polycorr(GParam *param, time_t time_start, time_t time_end)
{
  FILE *fp;
  int i;
  double diff_sec;

  fp = fopen(param->d_log_file, "w");
  fprintf(fp, "+--------------------------------------------+\n");
  fprintf(fp, "| Simulation details for yang_mills_polycorr |\n");
  fprintf(fp, "+--------------------------------------------+\n\n");

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

  fprintf(fp, "anisotropy: %d\n", param->d_anisotropic);
  fprintf(fp, "beta: %.10lf\n", param->d_beta);
  if (param->d_anisotropic != 0)
    fprintf(fp, "beta_t: %.10lf\n", param->d_beta_t);
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
  for (i = 0; i < NLEVELS; i++)
  {
    fprintf(fp, "%d ", param->d_ml_step[i]);
  }
  fprintf(fp, "\n");
  fprintf(fp, "updates for levels: ");
  for (i = 0; i < NLEVELS; i++)
  {
    fprintf(fp, "%d ", param->d_ml_upd[i]);
  }
  fprintf(fp, "\n");
  fprintf(fp, "dist_poly:  %d\n", param->d_dist_poly);
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
void print_parameters_t0(GParam *param, time_t time_start, time_t time_end)
{
  FILE *fp;
  int i;
  double diff_sec;

  fp = fopen(param->d_log_file, "w");
  fprintf(fp, "+--------------------------------------+\n");
  fprintf(fp, "| Simulation details for yang_mills_t0 |\n");
  fprintf(fp, "+--------------------------------------+\n\n");

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

  fprintf(fp, "gfstep:    %lf\n", param->d_gfstep);
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
void print_parameters_gf(GParam *param, time_t time_start, time_t time_end)
{
  FILE *fp;
  int i;
  double diff_sec;

  fp = fopen(param->d_log_file, "w");
  fprintf(fp, "+-------------------------------------------------+\n");
  fprintf(fp, "| Simulation details for yang_mills_gradient_flow |\n");
  fprintf(fp, "+-------------------------------------------------+\n\n");

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

  fprintf(fp, "gfstep:        %lf\n", param->d_gfstep);
  fprintf(fp, "num_gfsteps    %d\n", param->d_ngfsteps);
  fprintf(fp, "gf_meas_each   %d\n", param->d_gf_meas_each);
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

// print simulation parameters for the tracedef case
void print_parameters_tracedef(GParam const *const param, time_t time_start, time_t time_end, double acc)
{
  FILE *fp;
  int i;
  double diff_sec;

  fp = fopen(param->d_log_file, "w");
  fprintf(fp, "+--------------------------------------------+\n");
  fprintf(fp, "| Simulation details for yang_mills_tracedef |\n");
  fprintf(fp, "+--------------------------------------------+\n\n");

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

  fprintf(fp, "anisotropy: %d\n", param->d_anisotropic);
  fprintf(fp, "beta: %.10lf\n", param->d_beta);
  if (param->d_anisotropic != 0)
    fprintf(fp, "beta_t: %.10lf\n", param->d_beta_t);
#ifdef THETA_MODE
  fprintf(fp, "theta: %.10lf\n", param->d_theta);
#endif
  fprintf(fp, "h: %.10lf ", param->d_h[0]);
  for (i = 1; i < (int)floor(NCOLOR / 2.0); i++)
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
void print_parameters_tube_disc(GParam *param, time_t time_start, time_t time_end)
{
  FILE *fp;
  int i;
  double diff_sec;

  fp = fopen(param->d_log_file, "w");
  fprintf(fp, "+---------------------------------------------+\n");
  fprintf(fp, "| Simulation details for yang_mills_tube_disc |\n");
  fprintf(fp, "+---------------------------------------------+\n\n");

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

  fprintf(fp, "anisotropy: %d\n", param->d_anisotropic);
  fprintf(fp, "beta: %.10lf\n", param->d_beta);
  if (param->d_anisotropic != 0)
    fprintf(fp, "beta_t: %.10lf\n", param->d_beta_t);
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
  for (i = 0; i < NLEVELS; i++)
  {
    fprintf(fp, "%d ", param->d_ml_step[i]);
  }
  fprintf(fp, "\n");
  fprintf(fp, "updates for levels: ");
  for (i = 0; i < NLEVELS; i++)
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
void print_parameters_tube_conn(GParam *param, time_t time_start, time_t time_end)
{
  FILE *fp;
  int i;
  double diff_sec;

  fp = fopen(param->d_log_file, "w");
  fprintf(fp, "+---------------------------------------------+\n");
  fprintf(fp, "| Simulation details for yang_mills_tube_conn |\n");
  fprintf(fp, "+---------------------------------------------+\n\n");

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

  fprintf(fp, "anisotropy: %d\n", param->d_anisotropic);
  fprintf(fp, "beta: %.10lf\n", param->d_beta);
  if (param->d_anisotropic != 0)
    fprintf(fp, "beta_t: %.10lf\n", param->d_beta_t);
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
  for (i = 0; i < NLEVELS; i++)
  {
    fprintf(fp, "%d ", param->d_ml_step[i]);
  }
  fprintf(fp, "\n");
  fprintf(fp, "updates for levels: ");
  for (i = 0; i < NLEVELS; i++)
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
void print_parameters_tube_conn_long(GParam *param, time_t time_start, time_t time_end)
{
  FILE *fp;
  int i;
  double diff_sec;

  fp = fopen(param->d_log_file, "w");
  fprintf(fp, "+--------------------------------------------------+\n");
  fprintf(fp, "| Simulation details for yang_mills_tube_conn_long |\n");
  fprintf(fp, "+--------------------------------------------------+\n\n");

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

  fprintf(fp, "anisotropy: %d\n", param->d_anisotropic);
  fprintf(fp, "beta: %.10lf\n", param->d_beta);
  if (param->d_anisotropic != 0)
    fprintf(fp, "beta_t: %.10lf\n", param->d_beta_t);
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
  for (i = 0; i < NLEVELS; i++)
  {
    fprintf(fp, "%d ", param->d_ml_step[i]);
  }
  fprintf(fp, "\n");
  fprintf(fp, "updates for levels: ");
  for (i = 0; i < NLEVELS; i++)
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

#endif
