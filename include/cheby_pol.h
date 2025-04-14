#ifndef CHEBY_POLDEF
#define CHEBY_POLDEF

#include <stdio.h>
#include <stdlib.h>

#define MAX_CHEBY_DEG 10

typedef struct cheby_pol{
    int _max_deg;
    double coef[MAX_CHEBY_DEG];
}cheby_pol;

void init_cheby_pol(cheby_pol* p);
void cheby_pol_set_max_deg(cheby_pol* p, int deg);

void print_cheby_coef(cheby_pol* p, FILE* out);
void init_cheby_pol_from_file(cheby_pol* p, FILE* in);

double evaluate_cheby_pol(cheby_pol const * const p, double x);
double evaluete_cheby_pol_prime(cheby_pol const * const p, double x);
void cheby_pol_eval_and_prime(cheby_pol const * const p, double x, double* value, double* deriv);

double eval_cheby_pol_withgrad(cheby_pol const * const p, double const * const gradient);

void grad_coef_cheby_pol(double x, int len, double* gradient);
void grad_coef_cheby_pol_prime(double x, int len, double* gradient_prime);
void grad_coef_cheby_pol_and_prime(double x, int len, double* gradient, double* gradient_prime);

double x_to_negone_one(double val, double max, double min);
void project_to_fix_estrema(int len, double* gradient);

#endif