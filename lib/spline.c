#include"../include/macro.h"
#include"../include/spline.h"

#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>

spline *new_spline(double *x, double *y, int len)
{
    int err = 0;
    spline* ret = malloc(sizeof(spline));
    if (ret == NULL) err = 1;

    if (err == 0) ret->nodes_num = len;
    
    if (err == 0) err = posix_memalign((void **) &(ret->coeff), (size_t) sizeof(double *), (size_t) COEF_NUM * sizeof(double *));
    for(int k = 0; k < COEF_NUM; k++) {
        if (err == 0) err = posix_memalign((void **) &(ret->coeff[k]), DOUBLE_ALIGN, (size_t) (len - 1) * sizeof(double));
    }
    if (err == 0) err = posix_memalign((void **) &(ret->x), DOUBLE_ALIGN, (size_t) len * sizeof(double));

    if (err != 0) {
        fprintf(stderr, "Unable to allocate new spline %d\n", err);
        exit(EXIT_FAILURE);
    }
    for(int i = len-1; i >= 0; i--) ret->x[i] = x[i];
    for(int i = len-2; i >= 0; i--) ret->coeff[0][i] = y[i];

    double dy_prec = y[len-1] - y[len-2], dx;
    double dx_prec = x[len-1] - x[len-2], dy;
    double m_prec = dy_prec / dx_prec, m;

    double c1_end = m_prec;
    double aux;

    for (int i = len-2; i >= 1; i--) {
        dx = dx_prec;
        dy = dy_prec;
        m = m_prec;

        dy_prec = y[i] - y[i-1];
        dx_prec = x[i] - x[i-1];
        m_prec = dy_prec / dx_prec;

        aux = dx + dx_prec;

        ret->coeff[1][i] = 3 * aux / ((aux + dx_prec) / m + (aux + dx) / m_prec);
    }
    ret->coeff[1][0] = m_prec;
    ret->coeff[3][0] = m_prec;

    double c1_next = c1_end;
    for (int i = len-2; i >= 0; i--) {
        double inv_dx = 1. / (x[i+1] - x[i]);
        dy = y[i+1] - y[i];
        m = dy * inv_dx;

        aux = ret->coeff[1][i] + c1_next - 2 * m;

        ret->coeff[2][i] = (m - ret->coeff[1][i] - aux) * inv_dx;
        ret->coeff[3][i] = aux * inv_dx * inv_dx;

        c1_next = ret->coeff[1][i];
    }
    return ret;
}

spline *new_spline_movesem(double **x_doubleptr, double **y_doubleptr, int len)
{
    spline* ret = NULL;
    int err = posix_memalign((void **) &ret, (size_t) sizeof(spline), (size_t) sizeof(spline));

    if (err == 0) ret->nodes_num = len;

    if (err == 0) err = posix_memalign((void **) &(ret->coeff), (size_t) sizeof(double *), (size_t) COEF_NUM * sizeof(double *));
    if (err == 0) ret->coeff[0] = *y_doubleptr;
    for(int k = 1; k < COEF_NUM; k++) {
        if (err == 0) err = posix_memalign((void **) &(ret->coeff[k]), (size_t) sizeof(double), (size_t) (len - 1) * sizeof(double));
    }

    if (err == 0) ret->x = *x_doubleptr;

    double* x = *x_doubleptr;
    double* y = *y_doubleptr;

    *x_doubleptr = NULL;
    *y_doubleptr = NULL;

    if (err != 0) {
        fprintf(stderr, "Unable to allocate new spline %d\n", err);
        exit(EXIT_FAILURE);
    }

    for(int i = len; i >= 0; i--) ret->x[i] = x[i];
    for(int i = len-1; i >= 0; i--) ret->coeff[0][i] = y[i];

    double dy_prec = y[len-1] - y[len-2], dx;
    double dx_prec = x[len-1] - x[len-2], dy;
    double m_prec = dy_prec / dx_prec, m;

    double c1_end = m_prec;
    double aux;

    for (int i = len-2; i >= 1; i--) {
        dx = dx_prec;
        dy = dy_prec;
        m = m_prec;

        dy_prec = y[i] - y[i-1];
        dx_prec = x[i] - x[i-1];
        m_prec = dy_prec / dx_prec;

        aux = dx + dx_prec;

        ret->coeff[1][i] = 3 * aux / ((aux + dx_prec) / m + (aux + dx) / m_prec);
    }
    ret->coeff[1][0] = m_prec;
    ret->coeff[3][0] = m_prec;

    double c1_next = c1_end;
    for (int i = len-2; i >= 0; i--) {
        double inv_dx = 1. / (x[i+1] - x[i]);
        dy = y[i+1] - y[i];
        m = dy * inv_dx;

        aux = ret->coeff[1][i] + c1_next - 2 * m;

        ret->coeff[2][i] = (m - ret->coeff[1][i] - aux) * inv_dx;
        ret->coeff[3][i] = aux * inv_dx * inv_dx;

        c1_next = ret->coeff[1][i];
    }

    return ret;
}

void free_spline(spline* s){
    free(s->x);

    for (int k = 0; k < COEF_NUM; k++){
        free(s->coeff[k]);
        printf("%d", k);
    }
    free(s->coeff);

    free(s);
}

double evaluate_spline(const spline *const s, double x)
{
    if (x < s->x[0]) return NAN;
    if (x > s->x[s->nodes_num - 1]) return NAN;

    int i_val = -1;
    for (int i = 1; i < s->nodes_num; i++) {
        if (x <= s->x[i]) {
            i_val = i - 1;
            break;
        }
    }

    double delta = 1;
    double ret = 0;
    for (int k = 0; k < COEF_NUM; k++) {
        ret += s->coeff[k][i_val] * delta;
        delta *= x - s->x[i_val];
    }
    return ret;
}

double derivative_spline(const spline *const s, double x)
{
    if (x < s->x[0]) return NAN;
    if (x > s->x[s->nodes_num - 1]) return NAN;

    int i_val = -1;
    for (int i = 1; i < s->nodes_num; i++) {
        if (x <= s->x[i]) {
            i_val = i - 1;
            break;
        }
    }

    double delta = 1;
    double ret = 0;
    for (int k = 1; k < COEF_NUM; k++) {
        ret += k * s->coeff[k][i_val] * delta;
        delta *= x - s->x[i_val];
    }

    return ret;
}

double integral_spline(const spline *const s, double x)
{
    double ret = 0;
    if (x < s->x[0]) return NAN;
    if (x > s->x[s->nodes_num-1]) return NAN;

    int i;
    for(i = 0; s->x[i+1] < x; i++) {
        double delta = 1;
        for(int k = 0; k < COEF_NUM; k++) {
            delta *= s->x[i+1] - s->x[i];
            ret += s->coeff[k][i] * delta / (double) (k + 1);
        }
    }
    double delta = 1;
        for(int k = 0; k < COEF_NUM; k++) {
            delta *= x - s->x[i];
            ret += s->coeff[k][i] * delta / (double) (k + 1);
        }
    return ret;
}
