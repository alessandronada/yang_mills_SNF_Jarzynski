#include "../include/cheby_pol.h"

// initialize all coefficients to zero, but assuming max degree
void init_cheby_pol(cheby_pol *p)
{
    p->_max_deg = MAX_CHEBY_DEG;
    for (int k = 0; k < MAX_CHEBY_DEG; k++) p->coef[k] = 0;
}

// setting the degree to DEG
void cheby_pol_set_max_deg(cheby_pol *p, int deg)
{ if (deg > MAX_CHEBY_DEG) {
        fprintf(stderr, "Unable to set degree of Chebyshev polynomials to %d. Increase the maximum (MAX_CHEBY_DEG = %d) in cheby_pol.h",
                                                                          deg,                                      MAX_CHEBY_DEG);
        exit(EXIT_FAILURE);
    }
    p->_max_deg = deg;
}

void print_cheby_coef(cheby_pol *p, FILE *out)
{
    fprintf(out, "%d \n", p->_max_deg);
    for (int i = 0; i < p->_max_deg; i++) fprintf(out, "%lf ", p->coef[i]);
    fprintf(out, "\n");
}

void init_cheby_pol_from_file(cheby_pol *p, FILE *in)
{
    int read = fscanf(in, "%d \n", &(p->_max_deg));
    if (read != 1) {
        fprintf(stderr, "unable to read cheby_pol max deg \n");
        exit(EXIT_FAILURE);
    }
    if (p->_max_deg > MAX_CHEBY_DEG) {
        fprintf(stderr, "read cheby_pol max deg (%d) exceeds MAX_CHEBY_DEG (%d) \n", p->_max_deg, MAX_CHEBY_DEG);
        exit(EXIT_FAILURE);
    }
    for (int i = 0; i < p->_max_deg; i++) {
        if (read == 1) read = fscanf(in, "%lf ", &(p->coef[i]));
    }

    if (read != 1) {
        fprintf(stderr, "problems while reading cheby_pol coefficients \n");
        exit(EXIT_FAILURE);
    }

    for (int i = p->_max_deg; i < MAX_CHEBY_DEG; i++) {
        p->coef[i] = 0;
    }
}

// evaluate P in X
double evaluate_cheby_pol(cheby_pol const *const p, double x)
{   
    double ret = 0, T_now = 1, T_next = x, T_buff;
    for (int i = 0; i < p->_max_deg; i++) {
        ret += T_now * p->coef[i];

        T_buff = T_next;
        T_next = 2 * x * T_next - T_now;
        T_now = T_buff;
    }
    return ret;
}

// in practice: P->coef dot GRADIENT
double eval_cheby_pol_withgrad(cheby_pol const * const p, double const * const gradient)
{
    double ret = 0;
    for (int k = 0; k < p->_max_deg; k++) ret += gradient[k] * p->coef[k];
    return ret;
}

// The value of the base in X (T(X)), i.e. the gradient of p(X) wrt the parameters
void grad_coef_cheby_pol(double x, int len, double *gradient)
{
    double T_now = 1, T_next = x, T_buff;
    for (int i = 0; i < len; i++) {
        gradient[i] = T_now;

        T_buff = T_next;
        T_next = 2 * x * T_next - T_now;
        T_now = T_buff;
    }
}

// The value of the base and derivative of the base in X (T(X) and T'(X))
void grad_coef_cheby_pol_and_prime(double x, int len, double *gradient, double *gradient_prime)
/*
Computes the values and derivatives in X of the first LEN Chebyshv polynomials
and stores them in GRADIENT and GRADIENT_PRIME
(gradient with respect to parameters of the derivative of the linear
combination and its derivative)

GRADIENT and GRADIENT_PRIME are overwritten (can be uninitialized)
GRADIENT and GRADIENT_PRIME must point to at least LEN long arrays of double
*/
{
    double T_now = 1, U_prec = 0, T_next = x;

    gradient[0] = 1;
    gradient_prime[0] = 0;
    
    for (int i = 1; i < len; i++) {
        U_prec = x * U_prec + T_now;
        T_now = T_next;
        T_next = x * T_now - (1 - x * x) * U_prec;

        gradient[i] = T_now;
        gradient_prime[i] = i * U_prec;
    }
}

// Linear map from [MAX; MIN] into [-1; 1] evaluated in VAL
double x_to_negone_one(double val, double max, double min)
{
    return (2 * val - max - min) / (max - min);
}

// Project to the iperplane delta p(-1) = delta p(1) = 0, i.e. sum even = sum odd = 0
void project_to_fix_estrema(int len, double *gradient)
{
    double even_proj = 0, odd_proj = 0;
    for(int i = 1; i < len; i+=2) {
        even_proj += gradient[i - 1];
        odd_proj  += gradient[i];
    }
    if (len % 2) even_proj += gradient[len - 1];

    odd_proj  /= (double) (len / 2);
    even_proj /= (double) ((len + 1) / 2);
    for(int i = 1; i < len; i+=2) {
        gradient[i - 1] = gradient[i - 1] - even_proj;
        gradient[i] = gradient[i] - odd_proj;
    }
}
