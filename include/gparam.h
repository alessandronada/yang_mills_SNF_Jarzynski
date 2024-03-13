#ifndef GPARAM_H
#define GPARAM_H

#include"macro.h"

#include<stdio.h>
#include<stdlib.h>
#include<time.h>

typedef struct GParam {
  // lattice dimensions
  int d_size[STDIM];

  // simulation parameters
  double d_beta;
  double d_h[NCOLOR]; // parameters for the trace deformation
  double d_theta;
	
  // Jarzynski and parallel tempering parameters
  int d_defect_dir;				// defect boundary
  int d_L_defect[STDIM - 1];	// defect sizes

  // parallel tempering parameters
  int d_N_replica_pt;		// numbers of replica used in parallel tempering
  double* d_pt_bound_cond_coeff; // boundary conditions coefficients

  // Jarzynski parameters
  int d_J_evolutions; // number of out-of-equilibrium evolutions
  int d_J_relax;	// relax steps between evolutions
  int d_J_steps;	// steps in each evolution
	
	// hierarchical update (parallel tempering)
	int d_N_hierarc_levels; // number of hierarchical levels
	int *d_L_rect;          // d_L_rect is a vector of length d_N_hierarc_levels
	                        // d_L_rect[i] is the extension of the rectangle at the i-th hierarchical level
	int *d_N_sweep_rect;    // d_N_sweep_rect is vector of length d_N_hierarch_levels
												  // d_N_sweep_rect[i] is the number of sweep of the rectangle at the i-th hierarchical level

  // simulation details
  int d_sample;
  int d_thermal;
  int d_overrelax;
  int d_measevery;

  // initialization & saving
  int d_start;
  int d_saveconf_back_every;
  int d_saveconf_analysis_every;

  // for metropolis
  double d_epsilon_metro;

  // for cooling in measures
  int d_coolsteps;
  int d_coolrepeat;
	int d_chi_prime_meas;
	int d_topcharge_tprof_meas;

  // for gradient-flow evolution
  double d_gfstep;
	int d_ngfsteps;
	int d_gf_meas_each;

  // for multilevel
  int d_multihit;
  int d_ml_step[NLEVELS];
  int d_ml_upd[NLEVELS];
  int d_ml_level0_repeat;
  int d_dist_poly;
  int d_trasv_dist;
  int d_plaq_dir[2];

  // output file names
  char d_conf_file[STD_STRING_LENGTH];
  char d_data_file[STD_STRING_LENGTH];
	char d_chiprime_file[STD_STRING_LENGTH]; // print chi prime measures
	char d_topcharge_tprof_file[STD_STRING_LENGTH]; // print topological charge time correlator measures
  char d_log_file[STD_STRING_LENGTH];
  char d_ml_file[STD_STRING_LENGTH];
	char d_swap_acc_file[STD_STRING_LENGTH]; // print swap Metropolis acceptance
	char d_swap_tracking_file[STD_STRING_LENGTH];

  // random seed
  unsigned int d_randseed;

  // derived constants
  long d_volume;           // total volume
  double d_inv_vol;        // 1 / total volume
  long d_space_vol;        // spatial component of the volume
  double d_inv_space_vol;  // 1 / spatial component of the volume
	long d_volume_defect;		// volume of the defect (only for parallel tempering)
	int  d_n_grid;          // total grid points (only for multicanonic)
	
	// for multicanonic
	char d_topo_potential_file[STD_STRING_LENGTH];
	char d_multicanonic_acc_file[STD_STRING_LENGTH]; // print multicanonic Metropolis acceptance
	double d_grid_step;
	double d_grid_max;

} GParam;


void remove_white_line_and_comments(FILE *input);
void readinput(char *in_file, GParam *param);
void init_derived_constants(GParam *param);
void init_data_file(FILE **dataf, FILE **chiprimefilep, FILE **topchar_tprof_f, GParam const * const param);
void init_work_file(FILE **workfilep, GParam const * const param)
void free_hierarc_params(GParam *param);
void print_parameters_local(GParam const * const param, time_t time_start, time_t time_end);
void print_parameters_local_pt(GParam const * const param, time_t time_start, time_t time_end);
void print_parameters_local_jarzynski(GParam const * const param, time_t time_start, time_t time_end);
void print_parameters_local_pt_multicanonic(GParam const * const param, time_t time_start, time_t time_end);
void print_parameters_polycorr_long(GParam * param, time_t time_start, time_t time_end);
void print_parameters_polycorr(GParam * param, time_t time_start, time_t time_end);
void print_parameters_t0(GParam * param, time_t time_start, time_t time_end);
void print_parameters_gf(GParam * param, time_t time_start, time_t time_end);
void print_parameters_tracedef(GParam const * const param, time_t time_start, time_t time_end, double acc);
void print_parameters_tube_disc(GParam * param, time_t time_start, time_t time_end);
void print_parameters_tube_conn(GParam * param, time_t time_start, time_t time_end);
void print_parameters_tube_conn_long(GParam * param, time_t time_start, time_t time_end);

#endif
