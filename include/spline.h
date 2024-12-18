#ifndef MYSPLINEDEF
#define MYSPLINEDEF

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#define COEF_NUM 4

typedef struct spline{
    int nodes_num;
    double* x;
    double** coeff; // coeff[0][nodes_num], coeff[1:4][nodes_num-1]

    const char _free_flag;
} spline;

spline* new_spline(int len, double* xptr, double* yptr);
void init_spline(spline* s, int len, double* x, double* y);
void free_spline(spline* s);

void save_spline_onfile(const spline * const s, FILE* out);
void init_spline_fromfile(spline* s, FILE* in);

double evaluate_spline(const spline * const s, double x);
double derivative_spline(const spline * const s, double x);
double integral_spline(const spline * const s, double x);

#endif