#ifndef SUN_H
#define SUN_H

#include"macro.h"
#include"tens_prod.h"
#include"tens_prod_adj.h"

#include<complex.h>
#include<math.h>
#include<stdio.h>

typedef struct SuN {
   double complex comp[NCOLOR*NCOLOR] __attribute__((aligned(DOUBLE_ALIGN)));
} SuN;
//
//  the element [i][j] can be obtained by matrix.comp[m(i,j)] with m(i,j) defined in macro.h
//

typedef struct SuNAdj {
   #if NCOLOR!=1
     double comp[(NCOLOR*NCOLOR-1)*(NCOLOR*NCOLOR-1)] __attribute__((aligned(DOUBLE_ALIGN)));
   #else // this will never be used, is defined just to avoid warnings
     double comp[1] __attribute__((aligned(DOUBLE_ALIGN)));
   #endif

} SuNAdj;
//
//  the element [i][j] can be obtained by matrix.comp[madj(i,j)] with madj(i,j) defined in macro.h
//


// A=1
inline void one_SuN(SuN * restrict A)
  {
  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  #endif

  int i;

  for(i=0; i<NCOLOR*NCOLOR; i++)
     {
     A->comp[i]=0.0+0.0*I;
     }

  for(i=0; i<NCOLOR; i++)
     {
     A->comp[m(i,i)]=1.0+0.0*I;
     }
  }


// A=0
inline void zero_SuN(SuN * restrict A)
  {
  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  #endif

  int i;

  for(i=0; i<NCOLOR*NCOLOR; i++)
     {
     A->comp[i]=0.0+0.0*I;
     }
  }


// A=B
inline void equal_SuN(SuN * restrict A, SuN const * const restrict B)
  {
  #ifdef DEBUG
  if(A==B)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  __assume_aligned(&(B->comp), DOUBLE_ALIGN);
  #endif

  int i;

  for(i=0; i<NCOLOR*NCOLOR; i++)
     {
     A->comp[i]=B->comp[i];
     }
  }


// A=B^{dag}
inline void equal_dag_SuN(SuN * restrict A, SuN const * const restrict B)
  {
  #ifdef DEBUG
  if(A==B)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  __assume_aligned(&(B->comp), DOUBLE_ALIGN);
  #endif

  int i, j;

  for(i=0; i<NCOLOR; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        A->comp[m(i,j)]=conj(B->comp[m(j,i)]);
        }
     }
  }


// A+=B
inline void plus_equal_SuN(SuN * restrict A, SuN const * const restrict B)
  {
  #ifdef DEBUG
  if(A==B)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  __assume_aligned(&(B->comp), DOUBLE_ALIGN);
  #endif

  int i;

  for(i=0; i<NCOLOR*NCOLOR; i++)
     {
     A->comp[i]+=B->comp[i];
     }
  }


// A+=B^{dag}
inline void plus_equal_dag_SuN(SuN * restrict A, SuN const * const restrict B)
  {
  #ifdef DEBUG
  if(A==B)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  __assume_aligned(&(B->comp), DOUBLE_ALIGN);
  #endif

  int i, j;

  for(i=0; i<NCOLOR; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        A->comp[m(i,j)]+=conj(B->comp[m(j,i)]);
        }
     }
  }


// A-=B
inline void minus_equal_SuN(SuN * restrict A, SuN const * const restrict B)
  {
  #ifdef DEBUG
  if(A==B)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  __assume_aligned(&(B->comp), DOUBLE_ALIGN);
  #endif

  int i;

  for(i=0; i<NCOLOR*NCOLOR; i++)
     {
     A->comp[i]-=B->comp[i];
     }
  }


// A-=(r*B)
inline void minus_equal_times_real_SuN(SuN * restrict A, SuN const * const restrict B, double r)
  {
  #ifdef DEBUG
  if(A==B)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  __assume_aligned(&(B->comp), DOUBLE_ALIGN);
  #endif

  int i;

  for(i=0; i<NCOLOR*NCOLOR; i++)
     {
     A->comp[i]-=(r*B->comp[i]);
     }
  }


// A-=B^{dag}
inline void minus_equal_dag_SuN(SuN * restrict A, SuN const * const restrict B)
  {
  #ifdef DEBUG
  if(A==B)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  __assume_aligned(&(B->comp), DOUBLE_ALIGN);
  #endif

  int i, j;

  for(i=0; i<NCOLOR; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        A->comp[m(i,j)]-=conj(B->comp[m(j,i)]);
        }
     }
  }


// A=b*B+c*C
inline void lin_comb_SuN(SuN * restrict A,
                  double b, SuN const * const restrict B,
                  double c, SuN const * const restrict C)
  {
  #ifdef DEBUG
  if(A==B || A==C || B==C)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  __assume_aligned(&(B->comp), DOUBLE_ALIGN);
  __assume_aligned(&(C->comp), DOUBLE_ALIGN);
  #endif

  int i;

  for(i=0; i<NCOLOR*NCOLOR; i++)
     {
     A->comp[i]=b*(B->comp[i])+c*(C->comp[i]);
     }
  }


// A=b*B^{dag}+c*C
inline void lin_comb_dag1_SuN(SuN * restrict A,
                       double b, SuN const * const restrict B,
                       double c, SuN const * const restrict C)
  {
  #ifdef DEBUG
  if(A==B || A==C || B==C)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  __assume_aligned(&(B->comp), DOUBLE_ALIGN);
  __assume_aligned(&(C->comp), DOUBLE_ALIGN);
  #endif

  int i, j;

  for(i=0; i<NCOLOR; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        A->comp[m(i,j)]=b*conj(B->comp[m(j,i)])+c*(C->comp[m(i,j)]);
        }
     }
  }


// A=b*B+c*C^{dag}
inline void lin_comb_dag2_SuN(SuN * restrict A,
                       double b, SuN const * const restrict B,
                       double c, SuN const * const restrict C)
  {
  #ifdef DEBUG
  if(A==B || A==C || B==C)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  __assume_aligned(&(B->comp), DOUBLE_ALIGN);
  __assume_aligned(&(C->comp), DOUBLE_ALIGN);
  #endif

  int i, j;

  for(i=0; i<NCOLOR; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        A->comp[m(i,j)]=b*(B->comp[m(i,j)])+c*conj(C->comp[m(j,i)]);
        }
     }
  }


// A=b*B^{dag}+c*C^{dag}
inline void lin_comb_dag12_SuN(SuN * restrict A,
                        double b, SuN const * const restrict B,
                        double c, SuN const * const restrict C)
  {
  #ifdef DEBUG
  if(A==B || A==C || B==C)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  __assume_aligned(&(B->comp), DOUBLE_ALIGN);
  __assume_aligned(&(C->comp), DOUBLE_ALIGN);
  #endif

  int i, j;

  for(i=0; i<NCOLOR; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        A->comp[m(i,j)]=b*conj(B->comp[m(j,i)])+c*conj(C->comp[m(j,i)]);
        }
     }
  }


// A*=r
inline void times_equal_real_SuN(SuN * restrict A, double r)
  {
  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  #endif

  int i;

  for(i=0; i<NCOLOR*NCOLOR; i++)
     {
     A->comp[i]*=r;
     }
  }


// A*=r
inline void times_equal_complex_SuN(SuN * restrict A, double complex r)
  {
  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  #endif

  int i;

  for(i=0; i<NCOLOR*NCOLOR; i++)
     {
     A->comp[i]*=r;
     }
  }


// A*=B
inline void times_equal_SuN(SuN * restrict A, SuN const * const restrict B)
  {
  #ifdef DEBUG
  if(A==B)
   {
   fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
   exit(EXIT_FAILURE);
   }
  #endif

  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  __assume_aligned(&(B->comp), DOUBLE_ALIGN);
  #endif

  int i, j, k;
  double complex aux[NCOLOR] __attribute__((aligned(DOUBLE_ALIGN)));
  double complex sum;

  for(i=0; i<NCOLOR; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        aux[j]=A->comp[m(i,j)];
        }

     for(j=0; j<NCOLOR; j++)
        {
        sum=0.0+0.0*I;
        for(k=0; k<NCOLOR; k++)
           {
           sum+=aux[k]*(B->comp[m(k,j)]);
           }
        A->comp[m(i,j)]=sum;
        }
     }
  }


// A*=B^{dag}
inline void times_equal_dag_SuN(SuN * restrict A, SuN const * const restrict B)
  {
  #ifdef DEBUG
  if(A==B)
   {
   fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
   exit(EXIT_FAILURE);
   }
  #endif

  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  __assume_aligned(&(B->comp), DOUBLE_ALIGN);
  #endif

  int i, j, k;
  double complex aux[NCOLOR] __attribute__((aligned(DOUBLE_ALIGN)));
  double complex sum;

  for(i=0; i<NCOLOR; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        aux[j]=A->comp[m(i,j)];
        }

     for(j=0; j<NCOLOR; j++)
        {
        sum=0.0+0.0*I;
        for(k=0; k<NCOLOR; k++)
           {
           sum+=aux[k]*conj(B->comp[m(j,k)]);
           }
        A->comp[m(i,j)]=sum;
        }
     }
  }


// A=B*C
inline void times_SuN(SuN * restrict A,
                      SuN const * const restrict B,
                      SuN const * const restrict C)
  {
  #ifdef DEBUG
  if(A==B || A==C || B==C)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  __assume_aligned(&(B->comp), DOUBLE_ALIGN);
  __assume_aligned(&(C->comp), DOUBLE_ALIGN);
  #endif

  int i, j, k;
  double complex sum;

  for(i=0; i<NCOLOR; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        sum=0.0+0.0*I;
        for(k=0; k<NCOLOR; k++)
           {
           sum+=(B->comp[m(i,k)])*(C->comp[m(k,j)]);
           }
        A->comp[m(i,j)]=sum;
        }
     }
  }


// A=B^{dag}*C
inline void times_dag1_SuN(SuN * restrict A,
                    SuN const * const restrict B,
                    SuN const * const restrict C)
  {
  #ifdef DEBUG
  if(A==B || A==C || B==C)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  __assume_aligned(&(B->comp), DOUBLE_ALIGN);
  __assume_aligned(&(C->comp), DOUBLE_ALIGN);
  #endif

  int i, j, k;
  double complex sum;

  for(i=0; i<NCOLOR; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        sum=0.0+0.0*I;
        for(k=0; k<NCOLOR; k++)
           {
           sum+=conj(B->comp[m(k,i)])*(C->comp[m(k,j)]);
           }
        A->comp[m(i,j)]=sum;
        }
     }
  }


// A=B*C^{dag}
inline void times_dag2_SuN(SuN * restrict A,
                    SuN const * const restrict B,
                    SuN const * const restrict C)
  {
  #ifdef DEBUG
  if(A==B || A==C || B==C)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  __assume_aligned(&(B->comp), DOUBLE_ALIGN);
  __assume_aligned(&(C->comp), DOUBLE_ALIGN);
  #endif

  int i, j, k;
  double complex sum;

  for(i=0; i<NCOLOR; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        sum=0.0+0.0*I;
        for(k=0; k<NCOLOR; k++)
           {
           sum+=(B->comp[m(i,k)])*conj(C->comp[m(j,k)]);
           }
        A->comp[m(i,j)]=sum;
        }
     }
  }


// A=B^{dag}*C^{dag}
inline void times_dag12_SuN(SuN * restrict A,
                     SuN const * const restrict B,
                     SuN const * const restrict C)
  {
  #ifdef DEBUG
  if(A==B || A==C || B==C)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  __assume_aligned(&(B->comp), DOUBLE_ALIGN);
  __assume_aligned(&(C->comp), DOUBLE_ALIGN);
  #endif

  int i, j, k;
  double complex sum;

  for(i=0; i<NCOLOR; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        sum=0.0+0.0*I;
        for(k=0; k<NCOLOR; k++)
           {
           sum+=conj(B->comp[m(k,i)])*conj(C->comp[m(j,k)]);
           }
        A->comp[m(i,j)]=sum;
        }
     }
  }


// SU(N) random matrix
// generated a la Cabibbo Marinari with N(N-1)/2 SU(2) random matrices
void rand_matrix_SuN(SuN *A);


// generate a matrix in the algebra of SuN with gaussian
// random components in the base T_i such that Tr(T_iT_j)=delta_{ij}
void rand_algebra_gauss_matrix_SuN(SuN *A);


// l2 norm of the matrix
inline double norm_SuN(SuN const * const restrict A)
  {
  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  #endif

  int i;
  double aux, ris;

  ris=0.0;
  for(i=0; i<NCOLOR*NCOLOR; i++)
     {
     aux=cabs(A->comp[i]);
     ris+=aux*aux;
     }
  return sqrt(ris);
  }


// real part of the trace /N
inline double retr_SuN(SuN const * const restrict A)
  {
  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  #endif

  int i;
  double ris;
  double complex tr;

  tr=0.0+0.0*I;
  for(i=0; i<NCOLOR; i++)
     {
     tr+=A->comp[m(i,i)];
     }
  ris=creal(tr)/(double)NCOLOR;
  return ris;
  }


// imaginary part of the trace /N
inline double imtr_SuN(SuN const * const restrict A)
  {
  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  #endif

  int i;
  double ris;
  double complex tr;

  tr=0.0+0.0*I;
  for(i=0; i<NCOLOR; i++)
     {
     tr+=A->comp[m(i,i)];
     }
  ris=cimag(tr)/(double)NCOLOR;
  return ris;
  }


// LU decomposition with partial pivoting
void LU_SuN(SuN const * const A, SuN *ris, int *sign);


// determinant
inline complex double det_SuN(SuN const * const restrict A)
  {
  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  #endif

  #if NCOLOR==3
    complex double ris=0.0+0.0*I;

    ris+=(A->comp[m(0,0)])*(A->comp[m(1,1)])*(A->comp[m(2,2)]);
    ris+=(A->comp[m(1,0)])*(A->comp[m(2,1)])*(A->comp[m(0,2)]);
    ris+=(A->comp[m(2,0)])*(A->comp[m(0,1)])*(A->comp[m(1,2)]);
    ris-=(A->comp[m(2,0)])*(A->comp[m(1,1)])*(A->comp[m(0,2)]);
    ris-=(A->comp[m(1,0)])*(A->comp[m(0,1)])*(A->comp[m(2,2)]);
    ris-=(A->comp[m(0,0)])*(A->comp[m(2,1)])*(A->comp[m(1,2)]);

    return ris;
  #else
    int i;
    double complex ris;
    SuN lu;

    LU_SuN(A, &lu, &i);

    if(i>0)
      {
      ris=1.0+0.0*I;
      }
    else
     {
     ris=-1.0+0.0*I;
     }

    for(i=0; i<NCOLOR; i++)
       {
       ris*=(lu.comp[m(i,i)]);
       }

    return ris;
  #endif
  }


// gives 0 if the matrix is in SU(N) and 1 otherwise
int scheck_SuN(SuN const * const A);


// sunitarize
void unitarize_SuN(SuN *A);


// takes the traceless antihermitian part
inline void ta_SuN(SuN * restrict A)
  {
  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  #endif

  SuN aux, aux1;
  double complex trace;
  int i;

  equal_SuN(&aux, A);
  equal_dag_SuN(&aux1, A);
  minus_equal_SuN(&aux, &aux1);
  times_equal_real_SuN(&aux, 0.5); // now aux=(A-A^{dag})/2

  trace=aux.comp[m(0,0)];
  for(i=1; i<NCOLOR; i++)
     {
     trace+=aux.comp[m(i,i)];
     }
  trace/=(double)NCOLOR;

  for(i=0; i<NCOLOR; i++)
     {
     aux.comp[m(i,i)]-=trace;
     }

  equal_SuN(A, &aux);
  }

#if NCOLOR==3 // just to be sure
// https://arxiv.org/pdf/hep-lat/0311018 SEC III
inline void taexp_Su3(SuN * restrict A)
{
#ifdef __INTEL_COMPILER
   __assume_aligned(&(A->comp), DOUBLE_ALIGN);
#endif

   SuN aux, aux_sqr;
   equal_SuN(&aux, A);
   ta_SuN(&aux); // aux = 0.5 * (A - A^dagger - trace)
   times_equal_complex_SuN(&aux, -I); // aux is hermitian (eq. (2))

   double c0 = creal(det_SuN(&aux));
   int sign_c0;
   if (c0 >= 0) sign_c0 = 1;
   else {
      sign_c0 = -1;
      c0 = -c0;
   }

   equal_SuN(&aux_sqr, &aux);
   times_equal_SuN(&aux_sqr, &aux);

   // retr(.) = 1/3 Tr(.)
   double sqrt_c1_third = sqrt(0.5 * retr_SuN(&aux_sqr)); // sqrt(c1 / 3)
   double c0_max = 2. * pow(sqrt_c1_third, 3);
   double theta_third = acos(c0 / c0_max) / 3.;

   double u = sqrt_c1_third * cos(theta_third);
   double w = sqrt(3.) * sqrt_c1_third * sin(theta_third);
   double xi_0_w = fabs(w) > 0.05 ?
      sin(w) / w :
      1 - w * w / 6. * (1 - w * w / 20. * (1 - w * w / 42.));

   double h_real[3], h_imag[3];
   
   double cos_u = cos(u); // useful terms
   double cos_w = cos(w);
   double sin_u = sin(u);
   double cos_2u = cos(2 * u);
   double sin_2u = sin(2 * u);
   double u_sqr = u * u;
   double w_sqr = w * w;

   h_real[0] = (u_sqr - w_sqr) * cos_2u + (8 * u_sqr * cos_w) * cos_u + 2 * u * xi_0_w * (3 * u_sqr + w_sqr) * sin_u;
   h_imag[0] = (u_sqr - w_sqr) * sin_2u - (8 * u_sqr * cos_w) * sin_u + 2 * u * xi_0_w * (3 * u_sqr + w_sqr) * cos_u;

   h_real[1] = (2 * u) * cos_2u - (2 * u * cos_w) * cos_u + (3 * u_sqr - w_sqr) * xi_0_w * sin_u;
   h_imag[1] = (2 * u) * sin_2u + (2 * u * cos_w) * sin_u + (3 * u_sqr - w_sqr) * xi_0_w * cos_u;

   h_real[2] = cos_2u - cos_w * cos_u - 3 * u * xi_0_w * sin_u;
   h_imag[2] = sin_2u + cos_w * sin_u - 3 * u * xi_0_w * cos_u;

   // f_j(-c0) = (-1)^j f_j*(c0), eq(34)
   h_imag[0] *= sign_c0;
   h_real[1] *= sign_c0;
   h_imag[2] *= sign_c0;
   
   // h -> f;
   double denominator = 1. / (9 * u_sqr - w_sqr);
   for (int i = 0; i < 3; i++)
   {
      h_real[i] *= denominator;
      h_imag[i] *= denominator;
   }

   one_SuN(A);
   times_equal_complex_SuN(A, h_real[0] + h_imag[0] * I);

   times_equal_complex_SuN(&aux, h_real[1] + h_imag[1] * I);
   plus_equal_SuN(A, &aux);

   times_equal_complex_SuN(&aux_sqr, h_real[2] + h_imag[2] * I);
   plus_equal_SuN(A, &aux_sqr);
}

// Definitions from https://arxiv.org/pdf/2103.11965v3

// TP^i_j^k_l = A^k_j B^i_l
inline void oplus_SuN(TensProd * restrict TP, SuN const * const restrict A, SuN const *const restrict B) {
  #ifdef __INTEL_COMPILER
  __assume_aligned(&(B->comp), DOUBLE_ALIGN);
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  __assume_aligned(&(TP->comp), DOUBLE_ALIGN);
  #endif

  int i, j, k, l;

  for(i=0; i<NCOLOR; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        for(k=0; k<NCOLOR; k++)
           {
           for(l=0; l<NCOLOR; l++)
              {
              TP->comp[i][j][k][l]=A->comp[m(k,j)]*B->comp[m(i,l)];
              }
           }
        }
     }
}

// TP^i_j^k_l = A^i_j B^k_l
inline void otimes_SuN(TensProd * restrict TP, SuN const * const restrict A, SuN const *const restrict B) {
  #ifdef __INTEL_COMPILER
  __assume_aligned(&(B->comp), DOUBLE_ALIGN);
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  __assume_aligned(&(TP->comp), DOUBLE_ALIGN);
  #endif

  int i, j, k, l;

  for(i=0; i<NCOLOR; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        for(k=0; k<NCOLOR; k++)
           {
           for(l=0; l<NCOLOR; l++)
              {
              TP->comp[i][j][k][l]=A->comp[m(i,j)]*B->comp[m(k,l)];
              }
           }
        }
     }
}

inline void taexp_Su3_withderiv(SuN * restrict A, TensProd * restrict deriv)
{
#ifdef __INTEL_COMPILER
   __assume_aligned(&(A->comp), DOUBLE_ALIGN);
#endif

   SuN aux, aux_sqr;
   equal_SuN(&aux, A);
   ta_SuN(&aux); // aux = 0.5 * (A - A^dagger - trace)
   times_equal_complex_SuN(&aux, -I); // aux is hermitian (eq. (2))
   SuN Q; equal_SuN(&Q, &aux); // Q will not be changed
   SuN Q2; equal_SuN(&Q2, &aux); times_equal_SuN(&Q2, &aux); // Q^2 will not be changed

   double c0 = creal(det_SuN(&aux));
   int sign_c0;
   if (c0 >= 0) sign_c0 = 1;
   else {
      sign_c0 = -1;
      c0 = -c0;
   }

   equal_SuN(&aux_sqr, &Q2);
   //times_equal_SuN(&aux_sqr, &aux);

   // retr(.) = 1/3 Tr(.)
   double sqrt_c1_third = sqrt(0.5 * retr_SuN(&aux_sqr)); // sqrt(c1 / 3)
   double c0_max = 2. * pow(sqrt_c1_third, 3);
   double theta_third = acos(c0 / c0_max) / 3.;

   double u = sqrt_c1_third * cos(theta_third);
   double w = sqrt(3.) * sqrt_c1_third * sin(theta_third);
   double xi_0_w = fabs(w) > 0.05 ?
      sin(w) / w :
      1. - w * w / 6. * (1 - w * w / 20. * (1. - w * w / 42.));

   double h_real[3], h_imag[3];
   
   double cos_u = cos(u); // useful terms
   double cos_w = cos(w);
   double sin_u = sin(u);
   double cos_2u = cos(2. * u);
   double sin_2u = sin(2. * u);
   double u_sqr = u * u;
   double w_sqr = w * w;

   h_real[0] = (u_sqr - w_sqr) * cos_2u + (8. * u_sqr * cos_w) * cos_u + 2. * u * xi_0_w * (3. * u_sqr + w_sqr) * sin_u;
   h_imag[0] = (u_sqr - w_sqr) * sin_2u - (8. * u_sqr * cos_w) * sin_u + 2. * u * xi_0_w * (3. * u_sqr + w_sqr) * cos_u;

   h_real[1] = (2. * u) * cos_2u - (2. * u * cos_w) * cos_u + (3. * u_sqr - w_sqr) * xi_0_w * sin_u;
   h_imag[1] = (2. * u) * sin_2u + (2. * u * cos_w) * sin_u + (3. * u_sqr - w_sqr) * xi_0_w * cos_u;

   h_real[2] = cos_2u - cos_w * cos_u - 3. * u * xi_0_w * sin_u;
   h_imag[2] = sin_2u + cos_w * sin_u - 3. * u * xi_0_w * cos_u;
   
   // h -> f;
   double denominator = 1. / (9. * u_sqr - w_sqr);
   for (int i = 0; i < 3; i++)
   {
      h_real[i] *= denominator;
      h_imag[i] *= denominator;
   }

   // for the derivative
   double denominator2 = 0.5 * denominator * denominator;
   double xi_1_w = fabs(w) > 0.05 ?
      cos_w / w_sqr - sin(w) / (w * w_sqr) :
      -1./3. + w_sqr*(1./30. + w_sqr*(-1./840. + 1./45360.*w_sqr));

   double r_real[2][3], r_imag[2][3];
   //double **r_real = _r_real - 1; // never do r_real[0] !!!
   //double **r_imag = _r_imag - 1;

   double tmp_real, tmp_imag;

   // r_0^(1)
   tmp_real = (8. * u * cos_w + u * (3. * u_sqr + w_sqr) * xi_0_w);
   tmp_imag = (-4. * u_sqr * cos_w + (9. * u_sqr + w_sqr) * xi_0_w);
   r_real[0][0] = 2. * (u * cos_2u - (u_sqr - w_sqr) * sin_2u)
                + 2. * cos_u * tmp_real + 2. * sin_u * tmp_imag;
   r_imag[0][0] = 2. * (u * sin_2u + (u_sqr - w_sqr) * cos_2u)
                - 2. * sin_u * tmp_real + 2. * cos_u * tmp_imag;
   
   // r_1^(1)
   tmp_real = -2. * cos_w - (w_sqr - 3. * u_sqr) * xi_0_w;
   tmp_imag = 2. * u * cos_w + 6. * u * xi_0_w; 
   r_real[0][1] = 2. * (cos_2u - 2. * u * sin_2u) 
                + cos_u * tmp_real + sin_u * tmp_imag;
   r_imag[0][1] = 2. * (sin_2u + 2. * u * cos_2u)
                - sin_u * tmp_real + cos_u * tmp_imag;

   // r_2^(1)
   tmp_real = cos_w - 3. * xi_0_w;
   tmp_imag = 3. * u * xi_0_w;
   r_real[0][2] = -2. * sin_2u + sin_u * tmp_real - cos_u * tmp_imag;
   r_imag[0][2] = 2. * cos_2u + cos_u * tmp_real + sin_u * tmp_imag;

   // r_0^(2)
   tmp_real = cos_w + xi_0_w + 3. * u_sqr * xi_1_w;
   tmp_imag = 4. * u * xi_0_w;
   r_real[1][0] = -2. * cos_2u + 2. * u * sin_u * tmp_real - 2. * u * cos_u * tmp_imag;
   r_imag[1][0] = -2. * sin_2u + 2. * u * sin_u * tmp_imag + 2. * u * cos_u * tmp_real;

   // r_1^(2)
   tmp_real = cos_w + xi_0_w - 3. * u_sqr * xi_1_w;
   tmp_imag = 2. * u * xi_0_w;
   r_real[1][1] = -sin_u * tmp_real + cos_u * tmp_imag;
   r_imag[1][1] = -sin_u * tmp_imag - cos_u * tmp_real;

   // r_2^(2)
   tmp_real = xi_0_w;
   tmp_imag = -3. * u * xi_1_w;
   r_real[1][2] = cos_u * tmp_real + sin_u * tmp_imag;
   r_imag[1][2] = cos_u * tmp_imag - sin_u * tmp_real;

   double coeff1[3] = {-2. * (15. * u_sqr + w_sqr),
                   2. * u,
                   3. * u_sqr - w_sqr};
   double coeff2[3] = {-24. * u,
                   1.,
                   -3. * u};
   
   double b1_real[3], b1_imag[3];
   double b2_real[3], b2_imag[3];
   for (int j = 0; j < 3; j++){
      b1_real[j] = (coeff1[0] * h_real[j] + coeff1[1] * r_real[0][j] + coeff1[2] * r_real[1][j]) * denominator2;
      b1_imag[j] = (coeff1[0] * h_imag[j] + coeff1[1] * r_imag[0][j] + coeff1[2] * r_imag[1][j]) * denominator2;

      b2_real[j] = (coeff2[0] * h_real[j] + coeff2[1] * r_real[0][j] + coeff2[2] * r_real[1][j]) * denominator2;
      b2_imag[j] = (coeff2[0] * h_imag[j] + coeff2[1] * r_imag[0][j] + coeff2[2] * r_imag[1][j]) * denominator2;
   }

   // f_j(-c0) = (-1)^j f_j*(c0), eq(34)
   h_imag[0] *= sign_c0;
   h_real[1] *= sign_c0;
   h_imag[2] *= sign_c0;

   // b_ij(-c0) = (-1)^{i + j + 1} * b_ij*(c0)
   b1_imag[0] *= sign_c0;
   b1_real[1] *= sign_c0;
   b1_imag[2] *= sign_c0;
   b2_real[0] *= sign_c0;
   b2_imag[1] *= sign_c0;
   b2_real[2] *= sign_c0;

   one_SuN(A);
   times_equal_complex_SuN(A, h_real[0] + h_imag[0] * I);

   equal_SuN(&aux, &Q);
   times_equal_complex_SuN(&aux, h_real[1] + h_imag[1] * I);
   plus_equal_SuN(A, &aux);

   equal_SuN(&aux_sqr, &Q2);
   times_equal_complex_SuN(&aux_sqr, h_real[2] + h_imag[2] * I);
   plus_equal_SuN(A, &aux_sqr);
   
   SuN B1, B2;
   
   one_SuN(&B1); 
   times_equal_complex_SuN(&B1, b1_real[0] + I*b1_imag[0]); // B1 = b10

   equal_SuN(&aux, &Q);
   times_equal_complex_SuN(&aux, b1_real[1] + I*b1_imag[1]); 
   plus_equal_SuN(&B1, &aux); // B1 = b10 + b11 Q

   equal_SuN(&aux_sqr, &Q2);
   times_equal_complex_SuN(&aux_sqr, b1_real[2] + I*b1_imag[2]); 
   plus_equal_SuN(&B1, &aux_sqr);
   // now B1 = b10 + b11 Q + b12 Q^2

   //Fun! let's do it again
   one_SuN(&B2); 
   times_equal_complex_SuN(&B2, b2_real[0] + I*b2_imag[0]);

   equal_SuN(&aux, &Q);
   times_equal_complex_SuN(&aux, b2_real[1] + I*b2_imag[1]); 
   plus_equal_SuN(&B2, &aux);

   equal_SuN(&aux_sqr, &Q2);
   times_equal_complex_SuN(&aux_sqr, b2_real[2] + I*b2_imag[2]); 
   plus_equal_SuN(&B2, &aux_sqr);

   oplus_SuN(deriv, &Q, &B1); // deriv = Q oplus B1
   TensProd aux_TP; oplus_SuN(&aux_TP, &Q2, &B2); 
   plus_equal_TensProd(deriv, &aux_TP); // deriv = Q oplus B1 + Q^2 oplus B2
   
   // from now on aux is the identity
   one_SuN(&aux); 
   one_TensProd(&aux_TP);
   times_equal_complex_TensProd(&aux_TP, h_real[1] + I*h_imag[1]);
   plus_equal_TensProd(deriv, &aux_TP);
   
   otimes_SuN(&aux_TP, &Q, &aux);
   times_equal_complex_TensProd(&aux_TP, h_real[2] + I*h_imag[2]);
   plus_equal_TensProd(deriv, &aux_TP);

   otimes_SuN(&aux_TP, &aux, &Q); // can I merge this with previous?
   times_equal_complex_TensProd(&aux_TP, h_real[2] + I*h_imag[2]);
   plus_equal_TensProd(deriv, &aux_TP);
   // now deriv = Q oplus B1 + Q^2 oplus B2 + f1 * Id otimes Id + f2 * (Q otimes Id + Id otimes Q)
}
#endif // NCOLOR==3

// eponential of the traceless antihermitian part
inline void taexp_SuN(SuN * restrict A)
  {
  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  #endif

  SuN aux, uno, ris;

  equal_SuN(&aux, A);
  ta_SuN(&aux);

  one_SuN(&uno);

  // now aux is the traceless antihermitian part of the initial matrix
  // and we use
  // exp(x)=1+x(1+x/2(1+x/3*(1+x/4*(1+x/5*....

  equal_SuN(&ris, &aux);
  times_equal_real_SuN(&ris, 0.125);
  plus_equal_SuN(&ris, &uno);

  times_equal_SuN(&ris, &aux);
  times_equal_real_SuN(&ris, 0.142857142857142857142857);
  plus_equal_SuN(&ris, &uno);

  times_equal_SuN(&ris, &aux);
  times_equal_real_SuN(&ris, 0.16666666666666666666);
  plus_equal_SuN(&ris, &uno);

  times_equal_SuN(&ris, &aux);
  times_equal_real_SuN(&ris, 0.2);
  plus_equal_SuN(&ris, &uno);

  times_equal_SuN(&ris, &aux);
  times_equal_real_SuN(&ris, 0.25);
  plus_equal_SuN(&ris, &uno);

  times_equal_SuN(&ris, &aux);
  times_equal_real_SuN(&ris, 0.33333333333333333333);
  plus_equal_SuN(&ris, &uno);

  times_equal_SuN(&ris, &aux);
  times_equal_real_SuN(&ris, 0.5);
  plus_equal_SuN(&ris, &uno);

  times_equal_SuN(&ris, &aux);
  plus_equal_SuN(&ris, &uno);

  unitarize_SuN(&ris);
  equal_SuN(A, &ris);
  }



// return 0 if matrix is traceless antihermitian, 1 otherwise
inline int ta_check_SuN(SuN const * const restrict A)
  {
  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  #endif

  double complex aux;
  int i, j, ris;

  ris=0;

  aux=0.0;
  for(i=0; i<NCOLOR; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        aux+=(A->comp[m(i,j)]+conj(A->comp[m(j,i)]));
        }
     }
  if(cabs(aux)>MIN_VALUE) ris=1;

  aux=0.0;
  for(i=0; i<NCOLOR; i++)
     {
     aux+=A->comp[m(i,i)];
     }
  if(cabs(aux)>MIN_VALUE) ris=1;

  return ris;
  }


// exponential of a TA matrix
inline void exp_of_ta_SuN(SuN * restrict A)
  {
  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  #endif

  // we use
  // exp(x)=1+x(1+x/2(1+x/3*(1+x/4*(1+x/5*....

  #ifdef DEBUG
  if(ta_check_SuN(A)!=0)
    {
    fprintf(stderr, "Trying to exp. a non TA matrix! (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  SuN aux, uno;

  one_SuN(&uno);
  equal_SuN(&aux, A); // in aux the initial matrix is stored

  times_equal_real_SuN(A, 1.0/5.0);
  plus_equal_SuN(A, &uno);

  times_equal_SuN(A, &aux);
  times_equal_real_SuN(A, 1.0/4.0);
  plus_equal_SuN(A, &uno);

  times_equal_SuN(A, &aux);
  times_equal_real_SuN(A, 1.0/3.0);
  plus_equal_SuN(A, &uno);

  times_equal_SuN(A, &aux);
  times_equal_real_SuN(A, 1.0/2.0);
  plus_equal_SuN(A, &uno);

  times_equal_SuN(A, &aux);
  plus_equal_SuN(A, &uno);

  unitarize_SuN(A);
  }


// print on screen
void print_on_screen_SuN(SuN const * const A);


// print on file
void print_on_file_SuN(FILE *fp, SuN const * const A);


// print on binary file without changing endiannes
void print_on_binary_file_noswap_SuN(FILE *fp, SuN const * const A);


// print on binary file changing endiannes
void print_on_binary_file_swap_SuN(FILE *fp, SuN const * const A);


// print on binary file in bigendian
void print_on_binary_file_bigen_SuN(FILE *fp, SuN const * const A);


// read from file
void read_from_file_SuN(FILE *fp, SuN *A);


// read from binary file without changing endiannes
void read_from_binary_file_noswap_SuN(FILE *fp, SuN *A);


// read from binary file changing endianness
void read_from_binary_file_swap_SuN(FILE *fp, SuN *A);


// read from binary file written in bigendian
void read_from_binary_file_bigen_SuN(FILE *fp, SuN *A);


// initialize tensor product
inline void TensProd_init_SuN(TensProd * restrict TP, SuN const * const restrict A1, SuN const * const restrict A2)
  {
  #ifdef DEBUG
  if(A1==A2)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A1->comp), DOUBLE_ALIGN);
  __assume_aligned(&(A2->comp), DOUBLE_ALIGN);
  __assume_aligned(&(TP->comp), DOUBLE_ALIGN);
  #endif

  int i, j, k, l;

  for(i=0; i<NCOLOR; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        for(k=0; k<NCOLOR; k++)
           {
           for(l=0; l<NCOLOR; l++)
              {
              TP->comp[i][j][k][l]=conj(A1->comp[m(i,j)])*A2->comp[m(k,l)];
              }
           }
        }
     }
  }



// SuN B = (SuN A) star (TensProd T) or B^i_j = A^l_n T^l_j^i_n
inline void star_SuN_TensProd(SuN * restrict B, SuN const * const restrict A, TensProd const * const restrict TP)
{
#ifdef __INTEL_COMPILER
  __assume_aligned(&(B->comp), DOUBLE_ALIGN);
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  __assume_aligned(&(TP->comp), DOUBLE_ALIGN);
#endif
   int i, j, l, n;
   double complex sum;
   //zero_SuN(B);
   for (i = 0; i < NCOLOR; i++) {
      for (j = 0; j < NCOLOR; j++) {
         sum=0.0+0.0*I;
         for (l = 0; l < NCOLOR; l++) {
            for (n = 0; n < NCOLOR; n++) {
               sum += (A->comp[m(l, n)]) * (TP->comp[l][j][i][n]);
            }
         }
         B->comp[m(i, j)] = sum;
      }
   }
}

// TensProd A = (TensProd B) * (SuN M)
inline void times_rightSuN_TensProd(TensProd * restrict A,
                                    TensProd const * const restrict B,
                                    SuN const * const restrict M)
{
#ifdef DEBUG
   if (A == B) {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
   }
#endif

#ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  __assume_aligned(&(B->comp), DOUBLE_ALIGN);
  __assume_aligned(&(M->comp), DOUBLE_ALIGN);
#endif
   int i, j, k, l, n;   
   double complex sum;
   for (i = 0; i < NCOLOR; i++) {
      for (j = 0; j < NCOLOR; j++) {
         for (k = 0; k < NCOLOR; k++) {
            for (l = 0; l < NCOLOR; l++) {
               sum=0.0+0.0*I;
               for (n = 0; n < NCOLOR; n++) {
                  sum += (B->comp[i][j][k][n]) * (M->comp[m(n, l)]);
               }
               A->comp[i][j][k][l] = sum;
            }
         }
      }
   }
}

// TensProd A = (SuN M) * (TensProd B)
inline void times_leftSuN_TensProd(TensProd * restrict A,
                                    SuN const * const restrict M,
                                    TensProd const * const restrict B)
{
#ifdef DEBUG
   if (A == B) {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
   }
#endif
   int i, j, k, l, n;
   double complex sum;
   for (i = 0; i < NCOLOR; i++) {
      for (j = 0; j < NCOLOR; j++) {
         for (k = 0; k < NCOLOR; k++) {
            for (l = 0; l < NCOLOR; l++) {
               sum=0.0+0.0*I;
               for (n = 0; n < NCOLOR; n++) {
                  sum += (M->comp[m(i, n)]) * (B->comp[n][j][k][l]);
               }
               A->comp[i][j][k][l] = sum;
            }
         }
      }
   }
}

// convert the fundamental representation matrix B to the adjoint representation matrix A
inline void fund_to_adj_SuN(SuNAdj * restrict A, SuN const * const restrict B)
  {
  (void) A;
  (void) B;

  fprintf(stderr, "The function fund_to_adj_SuN still has to be written (%s, %d)\n", __FILE__, __LINE__);
  exit(EXIT_FAILURE);
  }


// initialize tensor product in the adjoint representation
// using two matrices in the fundamental representation
inline void TensProdAdj_init_SuN(TensProdAdj * restrict TP, SuN const * const restrict A1, SuN const * const restrict A2)
  {
  (void) TP;
  (void) A1;
  (void) A2;

  fprintf(stderr, "The function TensProd_adj_init_SuN still has to be written (%s, %d)\n", __FILE__, __LINE__);
  exit(EXIT_FAILURE);
  }


#endif
