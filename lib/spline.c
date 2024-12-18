#include"../include/macro.h"
#include"../include/spline.h"

#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>

// allocate new spline with LEN nodes. If not NULL, XPTR and YPTR are tight to nodes coordinate
spline *new_spline(int len, double* xptr, double* yptr)
{
    int err = 0;
    spline* ret = malloc(sizeof(spline));
    if (ret == NULL) err = 1;

    if (err == 0) ret->nodes_num = len;
    
    int free_flag = 0;
    if (xptr == NULL) {
        if (err == 0) err = posix_memalign((void **) &(ret->x), (size_t) sizeof(double), (size_t) len * sizeof(double));
        for (int i = 0; i < len; i++) ret->x[i] = 0;
    }
    else {
        ret->x = xptr;
        free_flag += 1;
    }
    
    if (err == 0) err = posix_memalign((void **) &(ret->coeff), (size_t) sizeof(double *), (size_t) COEF_NUM * sizeof(double *));
    if (yptr == NULL) {
        if (err == 0) err = posix_memalign((void **) &(ret->coeff[0]), (size_t) sizeof(double), (size_t) len * sizeof(double));
        if (err == 0) for (int i = 0; i < len; i++) ret->coeff[0][i] = 0;
    }
    else {
        ret->coeff[0] = yptr;
        free_flag += 2;
    }

    for(int k = 1; k < COEF_NUM; k++) {
        if (err == 0) err = posix_memalign((void **) &(ret->coeff[k]), (size_t) sizeof(double), (size_t) (len - 1) * sizeof(double));
        if (err == 0) for (int i = 0; i < len - 1; i++) ret->coeff[k][i] = 0;
    }

    if (err != 0) {
        fprintf(stderr, "Unable to allocate new spline %d\n", err);
        exit(EXIT_FAILURE);
    }

    * (char*) &(ret->_free_flag) = (char) free_flag;
    return ret;
}

// Initialize spline coefficients to inteprolate (X, Y) nodes (array of length LEN). If NULL X and Y are read from spline fields
void init_spline(spline* s, int len, double* x, double* y){
    if (s->nodes_num != len) {
        fprintf(stderr, "Len of points arrays does not match spline number of nodes! (%d != %d)\n", s->nodes_num, len);
        exit(EXIT_FAILURE);
    }

    if (x != NULL) for (int i = 0; i < len; i++) s->x[i] = x[i];
    else x = s->x;

    if (y != NULL) for (int i = 0; i < len; i++) s->coeff[0][i] = y[i];
    else y = s->coeff[0];

    double dy_prec = y[len-1] - y[len-2], dx;
    double dx_prec = x[len-1] - x[len-2], dy;
    double m_prec = dy_prec / dx_prec, m;

    double c1_next = m_prec; // saving the value for later
    double aux;

    for (int i = len-2; i >= 1; i--) {
        dx = dx_prec;
        dy = dy_prec;
        m = m_prec;

        dy_prec = y[i] - y[i-1];
        dx_prec = x[i] - x[i-1];
        m_prec = dy_prec / dx_prec;

        aux = dx + dx_prec;

        s->coeff[1][i] = 3 * aux / ((aux + dx_prec) / m + (aux + dx) / m_prec);
    }
    s->coeff[1][0] = m_prec;

    for (int i = len-2; i >= 0; i--) {
        double inv_dx = 1. / (x[i+1] - x[i]);
        dy = y[i+1] - y[i];
        m = dy * inv_dx;

        aux = s->coeff[1][i] + c1_next - 2 * m;

        s->coeff[2][i] = (m - s->coeff[1][i] - aux) * inv_dx;
        s->coeff[3][i] = aux * inv_dx * inv_dx;

        c1_next = s->coeff[1][i];
    }
}

// clean memory from S
void free_spline(spline* s){
    if (!(s->_free_flag & 1)) free(s->x);

    if (!(s->_free_flag & 2)) free(s->coeff[0]);
    for (int k = 1; k < COEF_NUM; k++){
        free(s->coeff[k]);
    }
    free(s->coeff);

    free(s);
}

// saving S on OUT as text
void save_spline_onfile(const spline *const s, FILE *out)
{
    fprintf(out, "%d\n", s->nodes_num);
    
    for (int i = 0; i < s->nodes_num; i++) fprintf(out, "%.12g ", s->x[i]);
    fprintf(out, "\n");
    
    for (int k = 0; k < COEF_NUM; k++) {
        for (int i = 0; i < s->nodes_num - 1; i++) fprintf(out, "%.12g ", s->coeff[k][i]);
        fprintf(out, "\n");
    }
}

// allocating and initializing a new spline reading coefficients from IN
void init_spline_fromfile(spline* s, FILE *in)
{
    int err = 0, len;

    if (err == 0) err = 1 - fscanf(in, "%d\n", &len);
    if (len != s->nodes_num) {
        fprintf(stderr, "Read number of nodes does not match spline's number of nodes (%d != %d)", len, s->nodes_num);
        exit(EXIT_FAILURE);
    }

    for(int i = 0; i < len; i++){
        if (err == 0) err = 1 - fscanf(in, "%lf ", &(s->x[i]));
    }
    for (int i = 0; i < len; i++) {
        if (err == 0) err = 1 - fscanf(in, "%lf ", &(s->coeff[0][i]));
    }
    for (int k = 1; k < COEF_NUM; k++) {
        for (int i = 0; i < len-1; i++) {
            if (err == 0) err = 1 - fscanf(in, "%lf ", &(s->coeff[k][i]));
        }
    }
    if (err != 0) {
        fprintf(stderr, "Error while reading spline drom file %d\n", err);
        exit(EXIT_FAILURE);
    }
}

// evaluate S on X
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

// evaluate S' (derivative of S) on X
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

// integrate S from S->x[0] (initial point) to X
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
