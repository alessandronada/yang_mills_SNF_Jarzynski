#ifndef MYSPLINEDEF
#define MYSPLINEDEF

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#define COEF_NUM 4

typedef struct spline{
    int nodes_num;
    double* x;
    double** coeff;
} spline;

spline* new_spline(double* x, double* y, int len);
spline* new_spline_movesem(double** x, double** y, int len);
void free_spline(spline* s);

double evaluate_spline(const spline * const s, double x);
double derivative_spline(const spline * const s, double x);
double integral_spline(const spline * const s, double x);

#endif