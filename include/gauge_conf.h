#ifndef GAUGE_CONF_H
#define GAUGE_CONF_H

#include"macro.h"

#include<complex.h>
#ifdef HASH_MODE
  #include<openssl/md5.h>
#endif
#include<stdio.h>

#include"gparam.h"
#include"geometry.h"
#include"su2.h"
#include"sun.h"
#include"tens_prod.h"
#include"tens_prod_adj.h"
#include"u1.h"


typedef struct Gauge_Conf {

  long update_index, evolution_index;

  GAUGE_GROUP **lattice;       // [volume] [STDIM]
  GAUGE_GROUP ***clover_array; // [volume] [STDIM] [STDIM]
	
	// for the parallel tempering of volume defect
	double **C;      // C [volume] [STDIM], this factor changes the boundary condition for the link on the defect
  int conf_label;  // save the label of the configuration to keep track of the swaps

  // for computing the polyakov loop correlator with multilevel
  TensProd ***ml_polycorr;   // [NLEVELS] [d_size[0]/d_ml_step[i]] [space_vol]
  GAUGE_GROUP **loc_poly;    // [d_size[0]/d_ml_step[NLEVELS-1]] [space_vol] auxilliary vector to be used in the multilevel

  // for computing the polyakov loop correlator in the adjoint rep. with multilevel
  TensProdAdj ***ml_polycorradj;   // [NLEVELS] [d_size[0]/d_ml_step[i]] [space_vol]

  // for the disconnected correlator for string width
  TensProd **ml_polyplaq;   // [NLEVELS] [only slice 0] [space_vol]
  double complex *loc_plaq;  // [only slice 0] [space_vol] auxilliary vector to be used in the multilevel

  // for the connected correlator for string width
  TensProd **ml_polyplaqconn;   // [NLEVELS] [only slice 0] [space_vol]
  GAUGE_GROUP *loc_plaqconn;    // [only slice 0][space_vol] auxilliary vector to be used in the multilevel
	
	// for multicanonic update (store running charge);
	double stored_topo_charge;
  } Gauge_Conf;
	
// to compute swap and multicanonic acceptances during parallel tempering evolution
typedef struct Acc_Utils {
  long *num_accepted_swap; // number of accepted swaps during parallel tempering
  long *num_swap;          // number of proposed swaps during parallel tempering
	
	long *num_accepted_metro_multicanonic; // number of accepted multicanonic Metropolis updates
	long *num_metro_multicanonic;          // number of proposed multicanonic Metropolis updates
  }	Acc_Utils;


// in gauge_conf_def.c
void init_gauge_conf_from_file_with_name(Gauge_Conf *GC,
                     GParam const * const param, char const * const filename);
void copy_gauge_conf_from_gauge_conf(Gauge_Conf *GC, Gauge_Conf const * const GC2, GParam const * const param);
void init_gauge_conf(Gauge_Conf *GC,
                     GParam const * const param);
void init_gauge_conf_replica(Gauge_Conf **GC,
															GParam const * const param);
void init_bound_cond(Gauge_Conf *GC,
										 GParam const * const param,
										 int const a);
void init_single_conf_bound_cond(Gauge_Conf *GC, GParam const * const param, double const Cval);
void set_bound_cond(Gauge_Conf *GC,	GParam const * const param,	double const Cval);
void free_replica(Gauge_Conf *GC,
									GParam const * const param);
void free_bound_cond(Gauge_Conf *GC,
										 GParam const * const param);
void read_gauge_conf_from_file_with_name(Gauge_Conf *GC,
                     GParam const * const param, char const * const filename);
void free_gauge_conf(Gauge_Conf *GC,
                    GParam const * const param);
void write_conf_on_file_with_name(Gauge_Conf const * const GC,
                                  GParam const * const param,
                                  char const * const namefile);
void write_conf_on_file(Gauge_Conf const * const GC,
                        GParam const * const param);
void write_conf_on_file_back(Gauge_Conf const * const GC,
                            GParam const * const param);
void write_evolution_conf_on_file(Gauge_Conf const* const GC, GParam const* const param, int obc);
void write_replica_on_file(Gauge_Conf const * const GC,
														GParam const * const param);
void write_replica_on_file_back(Gauge_Conf const * const GC,
																	GParam const * const param);
void init_gauge_conf_from_gauge_conf(Gauge_Conf *GC,
                                     Gauge_Conf const * const GC2,
                                     GParam const * const param);
void compute_md5sum_conf(char *res,        // the lenght is 2*MD5_DIGEST_LENGTH
                         Gauge_Conf const * const GC,
                         GParam const * const param);
void alloc_polycorr_stuff(Gauge_Conf *GC,
                          GParam const * const param);
void free_polycorr_stuff(Gauge_Conf *GC,
                         GParam const * const param);
void write_polycorr_on_file(Gauge_Conf const * const GC,
                            GParam const * const param,
                            int iteration);
void read_polycorr_from_file(Gauge_Conf const * const GC,
                             GParam const * const param,
                             int *iteration);
void compute_md5sum_polycorr(char *res,        // the lenght is 2*MD5_DIGEST_LENGTH
                             Gauge_Conf const * const GC,
                             GParam const * const param);
void alloc_polycorradj(Gauge_Conf *GC,
                       GParam const * const param);
void free_polycorradj(Gauge_Conf *GC,
                      GParam const * const param);
void alloc_tube_disc_stuff(Gauge_Conf *GC,
                           GParam const * const param);
void free_tube_disc_stuff(Gauge_Conf *GC,
                          GParam const * const param);
void alloc_tube_conn_stuff(Gauge_Conf *GC,
                           GParam const * const param);
void free_tube_conn_stuff(Gauge_Conf *GC,
                          GParam const * const param);
void write_tube_conn_stuff_on_file(Gauge_Conf const * const GC,
                                   GParam const * const param,
                                   int iteration);
void read_tube_conn_stuff_from_file(Gauge_Conf const * const GC,
                                    GParam const * const param,
                                    int *iteration);
void compute_md5sum_tube_conn_stuff(char *res, Gauge_Conf const * const GC, GParam const * const param);
void alloc_clover_array(Gauge_Conf *GC,
                        GParam const * const param);
void end_clover_array(Gauge_Conf *GC,
                      GParam const * const param);


// in gauge_conf_meas.c
double plaquettep(Gauge_Conf const * const GC,
                  Geometry const * const geo,
                  GParam const * const param,
                  long r,
                  int i,
                  int j);
double complex plaquettep_complex(Gauge_Conf const * const GC,
                                  Geometry const * const geo,
                                  GParam const * const param,
                                  long r,
                                  int i,
                                  int j);
void plaquettep_matrix(Gauge_Conf const * const GC,
                       Geometry const * const geo,
                       GParam const * const param,
                       long r,
                       int i,
                       int j,
                       GAUGE_GROUP *matrix);
void clover(Gauge_Conf const * const GC,
            Geometry const * const geo,
            GParam const * const param,
            long r,
            int i,
            int j,
            GAUGE_GROUP *M);
void plaquette(Gauge_Conf const * const GC,
               Geometry const * const geo,
               GParam const * const param,
               double *plaqs,
               double *plaqt);
void clover_disc_energy(Gauge_Conf const * const GC,
                        Geometry const * const geo,
                        GParam const * const param,
                        double *energy);
void polyakov(Gauge_Conf const * const GC,
              Geometry const * const geo,
              GParam const * const param,
              double *repoly,
              double *impoly);
void polyakov_adj(Gauge_Conf const * const GC,
                  Geometry const * const geo,
                  GParam const * const param,
                  double *repoly,
                  double *impoly);
void polyakov_with_tracedef(Gauge_Conf const * const GC,
                            Geometry const * const geo,
                            GParam const * const param,
                            double *repoly,
                            double *impoly);

void check_correlation_decay_cooling(Gauge_Conf const * const GC, Geometry const * const geo, GParam const * const param, double *ratio);
double sum_abs_topcharge_dens(Gauge_Conf const * const GC, Geometry const * const geo, GParam const * const param);

double loc_topcharge(Gauge_Conf const * const GC,
                     Geometry const * const geo,
                     GParam const * const param,
                     long r);
double topcharge(Gauge_Conf const * const GC,
                 Geometry const * const geo,
                 GParam const * const param);
void topcharge_timeslices(Gauge_Conf const * const GC,
                 Geometry const * const geo,
                 GParam const * const param, double *ris, int ncool, FILE*);
double topo_chi_prime(Gauge_Conf const * const GC,
                 Geometry const * const geo,
                 GParam const * const param);
void topo_obs_cooling(Gauge_Conf const * const GC,
                       Geometry const * const geo,
                       GParam const * const param,
                       double *charge,
											 double *chi_prime,
                       double *meanplaq, double *sum_q_timeslices, FILE *topchar_tprof_filep);
void topcharge_cooling(Gauge_Conf const * const GC,
                       Geometry const * const geo,
                       GParam const * const param,
                       double *charge,
                       double *meanplaq);
void loc_topcharge_corr(Gauge_Conf const * const GC,
                    Geometry const * const geo,
                    GParam const * const param,
                    int ncool,
                    int dist,
                    double *ris);
void perform_measures_localobs(Gauge_Conf *GC,
                               Geometry const * const geo,
                               GParam const * const param,
                               FILE *datafilep, FILE *chiprimefilep, FILE*);
void print_work(int index, double work, FILE *workfilep);
void perform_measures_localobs_with_tracedef(Gauge_Conf const * const GC,
                                             Geometry const * const geo,
                                             GParam const * const param,
                                             FILE *datafilep);
void optimize_multihit_polycorr(Gauge_Conf *GC,
                                Geometry const * const geo,
                                GParam const * const param,
                                FILE *datafilep);
void optimize_multilevel_polycorr(Gauge_Conf *GC,
                                  Geometry const * const geo,
                                  GParam const * const param,
                                  FILE *datafilep);
void perform_measures_polycorr(Gauge_Conf * GC,
                               Geometry const * const geo,
                               GParam const * const param,
                               FILE *datafilep);
void optimize_multihit_polycorradj(Gauge_Conf *GC,
                                   Geometry const * const geo,
                                   GParam const * const param,
                                   FILE *datafilep);
void optimize_multilevel_polycorradj(Gauge_Conf *GC,
                                     Geometry const * const geo,
                                     GParam const * const param,
                                     FILE *datafilep);
void perform_measures_polycorradj(Gauge_Conf *GC,
                                  Geometry const * const geo,
                                  GParam const * const param,
                                  FILE *datafilep);
void perform_measures_polycorr_long(Gauge_Conf * GC,
                                    GParam const * const param,
                                    FILE *datafilep);
void optimize_multilevel_polycorr_long(Gauge_Conf *GC,
                                       GParam const * const param,
                                       FILE *datafilep);
void perform_measures_polycorr_long(Gauge_Conf * GC,
                                    GParam const * const param,
                                    FILE *datafilep);
void perform_measures_tube_disc(Gauge_Conf *GC,
                                Geometry const * const geo,
                                GParam const * const param,
                                FILE *datafilep);
void perform_measures_tube_conn(Gauge_Conf *GC,
                                Geometry const * const geo,
                                GParam const * const param,
                                FILE *datafilep);
void perform_measures_tube_conn_long(Gauge_Conf *GC,
                                     GParam const * const param,
                                     FILE *datafilep);

// in gauge_conf_multilevel.c
void multihit(Gauge_Conf const * const GC,
              Geometry const * const geo,
              GParam const * const param,
              long r,
              int dir,
              int num_hit,
              GAUGE_GROUP *G);
void compute_local_poly(Gauge_Conf *GC,
                        Geometry const * const geo,
                        GParam const * const param);
void update_for_multilevel(Gauge_Conf * GC,
                           Geometry const * const geo,
                           GParam const * const param,
                           int level);
void multilevel_polycorr(Gauge_Conf *GC,
                          Geometry const * const geo,
                          GParam const * const param,
                          int dt);
void multilevel_polycorradj(Gauge_Conf * GC,
                            Geometry const * const geo,
                            GParam const * const param,
                            int dt);
void multilevel_polycorr_long(Gauge_Conf * GC,
                              Geometry const * const geo,
                              GParam const * const param,
                              int dt,
                              int iteration);
void compute_local_poly_and_plaq(Gauge_Conf *GC,
                                 Geometry const * const geo,
                                 GParam const * const param);
void multilevel_tube_disc(Gauge_Conf *GC,
                          Geometry const * const geo,
                          GParam const * const param,
                          int dt);
void compute_local_poly_plaq_and_plaqconn(Gauge_Conf *GC,
                                          Geometry const * const geo,
                                          GParam const * const param);
void multilevel_tube_conn(Gauge_Conf * GC,
                          Geometry const * const geo,
                          GParam const * const param,
                          int dt);
void multilevel_tube_conn_long(Gauge_Conf * GC,
                               Geometry const * const geo,
                               GParam const * const param,
                               int dt,
                               int iteration);

// in gauge_conf_upd.c
void calcstaples_wilson(Gauge_Conf const * const GC,
                        Geometry const * const geo,
                        GParam const * const gparam,
                        long r,
                        int i,
                        GAUGE_GROUP *M);
void calcstaples_tracedef(Gauge_Conf const * const GC,
                          Geometry const * const geo,
                          GParam const * const param,
                          long r,
                          int i,
                          GAUGE_GROUP * M);
void compute_clovers(Gauge_Conf const * const GC,
                     Geometry const * const geo,
                     GParam const * const param,
                     int direction);
void calcstaples_with_topo(Gauge_Conf const * const GC,
                           Geometry const * const geo,
                           GParam const * const param,
                           long r,
                           int i,
                           GAUGE_GROUP *M);
void compute_clovers_replica(Gauge_Conf const * const GC, Geometry const * const geo,
													GParam const * const param, int dir);
void compute_clovers_replica_rect(Gauge_Conf const * const GC, Geometry const * const geo,
													GParam const * const param, int dir, Rectangle const * const clover_rectangle);
void calcstaples_wilson_with_defect(Gauge_Conf const * const GC,
                        Geometry const * const geo,
                        GParam const * const param,
                        long r,
                        int i,
                        GAUGE_GROUP *M);
void calcstaples_with_topo_with_defect(Gauge_Conf const * const GC,
                           Geometry const * const geo,
                           GParam const * const param,
                           long r,
                           int i,
                           GAUGE_GROUP *M);
void heatbath(Gauge_Conf * GC,
              Geometry const * const geo,
              GParam const * const param,
              long r,
              int i);
void overrelaxation(Gauge_Conf * GC,
                    Geometry const * const geo,
                    GParam const * const param,
                    long r,
                    int i);
int metropolis(Gauge_Conf *GC,
               Geometry const * const geo,
               GParam const * const param,
               long r,
               int i);
int metropolis_with_tracedef(Gauge_Conf *GC,
                             Geometry const * const geo,
                             GParam const * const param,
                             long r,
                             int i);
void heatbath_with_defect(Gauge_Conf *GC,
              Geometry const * const geo,
              GParam const * const param,
              long r,
              int i);
void overrelaxation_with_defect(Gauge_Conf *GC,
                    Geometry const * const geo,
                    GParam const * const param,
                    long r,
                    int i);
void update(Gauge_Conf *GC,
            Geometry const * const geo,
            GParam const * const param);
void update_with_defect(Gauge_Conf * GC,
												Geometry const * const geo,
												GParam const * const param);
void update_rectangle_with_defect(Gauge_Conf *GC, Geometry const * const geo, GParam const * const param,
											Rectangle const * const most_update, Rectangle const * const clover_rectangle);
void single_conf_hierarchical_update_rectangle_with_defect(Gauge_Conf *GC, Geometry const * const geo, GParam const * const param,
	int const hierarc_level, Rectangle const * const most_update,
	Rectangle const * const clover_rectangle);
void hierarchical_update_rectangle_with_defect(Gauge_Conf *GC, Geometry const * const geo, GParam const * const param,
															int const hierarc_level, Rectangle const * const most_update,
															Rectangle const * const clover_rectangle,
															Rectangle const * const swap_rectangle, Acc_Utils *acc_counters);
void single_conf_hierarchical_update(Gauge_Conf *GC, Geometry const * const geo, GParam const * const param,
	Rectangle const * const most_update, Rectangle const * const clover_rectangle);
void parallel_tempering_with_hierarchical_update(Gauge_Conf *GC, Geometry const * const geo, GParam const * const param,
												 Rectangle const * const most_update, Rectangle const * const clover_rectangle,
												 Rectangle const * const swap_rectangle, Acc_Utils *acc_counters);
void update_with_trace_def(Gauge_Conf *GC,
                           Geometry const * const geo,
                           GParam const * const param,
                           double *acc);
void cooling(Gauge_Conf *GC,
             Geometry const * const geo,
             GParam const * const param,
             int n);
void gradflow_RKstep(Gauge_Conf *GC,
                     Gauge_Conf *helper1,
                     Gauge_Conf *helper2,
                     Geometry const * const geo,
                     GParam const *const param,
                     double dt);
void ape_smearing(Gauge_Conf *GC,
                  Geometry const * const geo,
                  GParam const *const param,
                  double alpha,
                  int n);
									
// in gauge_conf_paral_temp.c
void swap(Gauge_Conf *GC, Geometry const * const geo, GParam const * const param,
				 Rectangle const * const swap_rectangle, Acc_Utils *acc_counters);
double delta_action_swap(Gauge_Conf const * const GC, Geometry const * const geo, GParam const * const param,
                         long const r, int const i, int const j, int const a, int const b);
double compute_defect_action(Gauge_Conf const * const GC, Geometry const * const geo, GParam const * const param);
void metropolis_single_swap(Gauge_Conf *GC, int const a, int const b, double const p, Acc_Utils *acc_counters);
void conf_translation(Gauge_Conf *GC, Geometry const * const geo, GParam const * const param);	
void init_swap_acc_arrays(Acc_Utils *acc_counters, GParam const * const param);
void end_swap_acc_arrays(Acc_Utils *acc_counters, GParam const * const param);
void print_acceptances(Acc_Utils const * const acc_counters, GParam const * const param);
void init_swap_track_file(FILE **swaptrackfilep, GParam const * const param);
void print_conf_labels(FILE *fp, Gauge_Conf const * const GC, GParam const * const param);

// in gauge_conf_multicanonic.c
void read_topo_potential(double **, GParam const * const);
double compute_topo_potential(double const, double const * const, GParam const * const);
void init_multicanonic_acc_file(FILE **, GParam const * const);
void init_multicanonic_acc_arrays(Acc_Utils *, GParam const * const);
void end_multicanonic_acc_arrays(Acc_Utils *);
void print_multicanonic_acceptance(Gauge_Conf const * const, GParam const * const, Acc_Utils const * const, FILE *);
void multicanonic_parallel_tempering_with_hierarchical_update(Gauge_Conf *, Geometry const * const, GParam const * const, Rectangle const * const, Rectangle const * const,
																															Rectangle const * const, Acc_Utils *, double const * const, FILE *);
void multicanonic_update_with_defect(Gauge_Conf *, Geometry const * const, GParam const * const, double const * const, Acc_Utils *acc_counters);
void multicanonic_update_rectangle_with_defect(Gauge_Conf *, Geometry const * const, GParam const * const,
																								Rectangle const * const, Rectangle const * const, double const * const, Acc_Utils *);
void multicanonic_hierarchical_update_rectangle_with_defect(Gauge_Conf *, Geometry const * const, GParam const * const,
														 int const, Rectangle const * const, Rectangle const * const,
														 Rectangle const * const, Acc_Utils *, double const * const);
int multicanonic_overrelaxation_with_defect(Gauge_Conf *, Geometry const * const, GParam const * const,
																							long, int, double const * const);
int multicanonic_heatbath_with_defect(Gauge_Conf *, Geometry const * const, GParam const * const,
																				long, int, double const * const);
void init_topo_charge(Gauge_Conf *, Geometry const * const, GParam const * const);
void refresh_topo_charge_replica(Gauge_Conf *, Geometry const * const, GParam const * const);
void refresh_topo_charge(Gauge_Conf *, Geometry const * const, GParam const * const);
void compute_single_clover_insertion(GAUGE_GROUP *, Gauge_Conf const * const, Geometry const * const, GParam const * const param, long, int, int);
void compute_topostaple_alone(Gauge_Conf const * const, Geometry const * const, GParam const * const, long, int, GAUGE_GROUP *);								 
double delta_Q_upd(Gauge_Conf const * const, Geometry const * const, GParam const * const, long, int, GAUGE_GROUP);
double metropolis_prob_multicanonic(double const, double const, GParam const * const, double const * const);
int multicanonic_Metropolis_step_single_link(Gauge_Conf *, Geometry const * const, GParam const * const, long, int, double const * const, GAUGE_GROUP);

#endif
