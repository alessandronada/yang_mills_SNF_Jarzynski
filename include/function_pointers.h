#ifndef FUNCTION_POINTERS_H
#define FUNCTION_POINTERS_H

#include"macro.h"
#include"su2.h"
#include"su2_upd.h"
#include"sun.h"
#include"sun_upd.h"
#include"tens_prod.h"
#include"tens_prod_adj.h"
#include"u1.h"
#include"u1_upd.h"

#include<complex.h>
#include<stdio.h>

// declarations only: the pointers are defined (and initialized) in lib/function_pointers.c;
// without extern every translation unit would define them, which gcc >= 10 rejects at link time

extern void (*one)(GAUGE_GROUP *A);         // A=1
extern void (*zero)(GAUGE_GROUP *A);        // A=0

extern void (*equal)(GAUGE_GROUP *A,
                     GAUGE_GROUP const * const B);      // A=B
extern void (*equal_dag)(GAUGE_GROUP *A,
                         GAUGE_GROUP const * const B);  // A=B^{dag}

extern void (*plus_equal)(GAUGE_GROUP *A,
                          GAUGE_GROUP const * const B);     // A+=B
extern void (*plus_equal_dag)(GAUGE_GROUP *A,
                              GAUGE_GROUP const * const B); // A+=B^{dag}

extern void (*minus_equal)(GAUGE_GROUP *A,
                           GAUGE_GROUP const * const B);     // A-=B
extern void (*minus_equal_times_real)(GAUGE_GROUP *A,
                                      GAUGE_GROUP const * const B, double r);   // A-=(r*B)
extern void (*minus_equal_dag)(GAUGE_GROUP *A,
                               GAUGE_GROUP const * const B); // A-=B^{dag}

extern void (*lin_comb)(GAUGE_GROUP *A,
                        double b, GAUGE_GROUP const * const B,
                        double c, GAUGE_GROUP const * const C);       // A=b*B+c*C
extern void (*lin_comb_dag1)(GAUGE_GROUP *A,
                             double b, GAUGE_GROUP const * const B,
                             double c, GAUGE_GROUP const * const C);  // A=b*B^{dag}+c*C
extern void (*lin_comb_dag2)(GAUGE_GROUP *A,
                             double b, GAUGE_GROUP const * const B,
                             double c, GAUGE_GROUP const * const C);  // A=b*B+c*C^{dag}
extern void (*lin_comb_dag12)(GAUGE_GROUP *A,
                              double b, GAUGE_GROUP const * const B,
                              double c, GAUGE_GROUP const * const C); // A=b*B^{dag}+c*C^{dag}

extern void (*times_equal_real)(GAUGE_GROUP *A, double r); // A*=r
extern void (*times_equal_complex)(GAUGE_GROUP *A, double complex r); // A*=r

extern void (*times_equal)(GAUGE_GROUP *A,
                           GAUGE_GROUP const * const B);     // A*=B
extern void (*times_equal_dag)(GAUGE_GROUP *A,
                               GAUGE_GROUP const *B); // A*=B^{dag}

extern void (*times)(GAUGE_GROUP *A,
                     GAUGE_GROUP const * const B,
                     GAUGE_GROUP const * const C);       // A=B*C
extern void (*times_dag1)(GAUGE_GROUP *A,
                          GAUGE_GROUP const * const B,
                          GAUGE_GROUP const * const C);  // A=B^{dag}*C
extern void (*times_dag2)(GAUGE_GROUP *A,
                          GAUGE_GROUP const * const B,
                          GAUGE_GROUP const * const C);  // A=B*C^{dag}
extern void (*times_dag12)(GAUGE_GROUP *A,
                           GAUGE_GROUP const * const B,
                           GAUGE_GROUP const * const C); // A=B^{dag}*C^{dag}

extern void (*rand_matrix)(GAUGE_GROUP *A);

extern double (*norm)(GAUGE_GROUP const * const A);

extern double (*retr)(GAUGE_GROUP const * const A);
extern double (*imtr)(GAUGE_GROUP const * const A);

extern void (*unitarize)(GAUGE_GROUP *A);
extern void (*ta)(GAUGE_GROUP *A);
extern void (*taexp)(GAUGE_GROUP *A);

extern void (*print_on_screen)(GAUGE_GROUP const * const A);
extern void (*print_on_file)(FILE *fp, GAUGE_GROUP const * const A);
extern void (*print_on_binary_file_bigen)(FILE *fp, GAUGE_GROUP const * const A);
extern void (*read_from_file)(FILE *fp, GAUGE_GROUP *A);
extern void (*read_from_binary_file_bigen)(FILE *fp, GAUGE_GROUP *A);

extern void (*TensProd_init)(TensProd *TP, GAUGE_GROUP const * const A1, GAUGE_GROUP const * const A2);
extern void (*TensProdAdj_init)(TensProdAdj *TP, GAUGE_GROUP const * const A1, GAUGE_GROUP const * const A2);

extern void (*single_heatbath)(GAUGE_GROUP *link, GAUGE_GROUP const * const staple, GParam const * const param);
extern void (*single_overrelaxation)(GAUGE_GROUP *link, GAUGE_GROUP const * const staple);
extern void (*cool)(GAUGE_GROUP *link, GAUGE_GROUP const * const staple);

#endif
