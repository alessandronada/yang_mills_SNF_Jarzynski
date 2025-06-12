#ifndef TENS_PROD_H
#define TENS_PROD_H

#include"macro.h"

#include<complex.h>
#include<stdio.h>
#include<stdlib.h>

// see Luscher Weisz JHEP 0109 p. 010 (2001)   (hep-lat/0108014)

typedef struct TensProd {
   double complex comp[NCOLOR][NCOLOR][NCOLOR][NCOLOR] __attribute__((aligned(DOUBLE_ALIGN)));
} TensProd;


// initialize to zero
inline void zero_TensProd(TensProd * restrict A)
  {
  int i0, i1, i2, i3;

  for(i0=0; i0<NCOLOR; i0++)
     {
     for(i1=0; i1<NCOLOR; i1++)
        {
        for(i2=0; i2<NCOLOR; i2++)
           {
           for(i3=0; i3<NCOLOR; i3++)
              {
              A->comp[i0][i1][i2][i3]=0.0+0.0*I;
              }
           }
        }
     }
  }


// initialize to one
inline void one_TensProd(TensProd * restrict A)
  {
  int i0, i1, i2, i3;

  for(i0=0; i0<NCOLOR; i0++)
     {
     for(i1=0; i1<NCOLOR; i1++)
        {
        for(i2=0; i2<NCOLOR; i2++)
           {
           for(i3=0; i3<NCOLOR; i3++)
              {
              A->comp[i0][i1][i2][i3]=0.0+0.0*I;
              }
           }
        }
     }

  for(i0=0; i0<NCOLOR; i0++)
     {
     for(i1=0; i1<NCOLOR; i1++)
        {
        A->comp[i0][i0][i1][i1]=1.0+0.0*I;
        }
     }
  }


// A=B
inline void equal_TensProd(TensProd * restrict A, TensProd const * const restrict B)
  {
  #ifdef DEBUG
  if(A==B)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  int i0, i1, i2, i3;

  for(i0=0; i0<NCOLOR; i0++)
     {
     for(i1=0; i1<NCOLOR; i1++)
        {
        for(i2=0; i2<NCOLOR; i2++)
           {
           for(i3=0; i3<NCOLOR; i3++)
              {
              A->comp[i0][i1][i2][i3]=B->comp[i0][i1][i2][i3];
              }
           }
        }
     }
  }


// A*=r
inline void times_equal_real_TensProd(TensProd * restrict A, double r)
  {
  int i0, i1, i2, i3;

  for(i0=0; i0<NCOLOR; i0++)
     {
     for(i1=0; i1<NCOLOR; i1++)
        {
        for(i2=0; i2<NCOLOR; i2++)
           {
           for(i3=0; i3<NCOLOR; i3++)
              {
              A->comp[i0][i1][i2][i3]*=r;
              }
           }
        }
     }
  }


// A*=r
inline void times_equal_complex_TensProd(TensProd * restrict A, double complex r)
  {
  int i0, i1, i2, i3;

  for(i0=0; i0<NCOLOR; i0++)
     {
     for(i1=0; i1<NCOLOR; i1++)
        {
        for(i2=0; i2<NCOLOR; i2++)
           {
           for(i3=0; i3<NCOLOR; i3++)
              {
              A->comp[i0][i1][i2][i3]*=r;
              }
           }
        }
     }
  }


// A+=B
inline void plus_equal_TensProd(TensProd * restrict A, TensProd const * const restrict B)
  {
  #ifdef DEBUG
  if(A==B)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  int i0, i1, i2, i3;

  for(i0=0; i0<NCOLOR; i0++)
     {
     for(i1=0; i1<NCOLOR; i1++)
        {
        for(i2=0; i2<NCOLOR; i2++)
           {
           for(i3=0; i3<NCOLOR; i3++)
              {
              A->comp[i0][i1][i2][i3]+=B->comp[i0][i1][i2][i3];
              }
           }
        }
     }
  }


// A=B*C
inline void times_TensProd(TensProd * restrict A,
                           TensProd const * const restrict B,
                           TensProd const * const restrict C)
  {
  #ifdef DEBUG
  if(A==B || A==C || B==C)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  int i, j, k, l;
  int a,b;

  double complex sum;

  for(i=0; i<NCOLOR; i++)
     {
     for(j=0; j<NCOLOR; j++)
        {
        for(k=0; k<NCOLOR; k++)
           {
           for(l=0; l<NCOLOR; l++)
              {
              sum=0.0+0.0*I;
              for(a=0; a<NCOLOR; a++)
                 {
                 for(b=0; b<NCOLOR; b++)
                    {
                    sum+=B->comp[i][a][k][b] * C->comp[a][j][b][l];
                    }
                 }
              A->comp[i][j][k][l]=sum;
              }
           }
        }
     }
  }


// A*=B
inline void times_equal_TensProd(TensProd * restrict A, TensProd const * const restrict B)
  {
  #ifdef DEBUG
  if(A==B )
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  TensProd tmp __attribute__((aligned(DOUBLE_ALIGN)));

  equal_TensProd(&tmp, A);
  times_TensProd(A, &tmp, B);
  }

// in https://arxiv.org/pdf/2103.11965v3 A = B star C
// A^i_j^k_l = B^i_n^m_l C^n_j^k_m
inline void star_TensProd(TensProd * restrict A,
                          TensProd const * const restrict B,
                          TensProd const * const restrict C)
{
#ifdef DEBUG
   if(A==B || A==C || B==C)
      {
      fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
      exit(EXIT_FAILURE);
      }
#endif
   // zero_TensProd(A);
   int i, j, k, l, m, n;
   double complex sum;
   
   // this is not pretty
   for (i = 0; i < NCOLOR; i++) {
      for (j = 0; j < NCOLOR; j++) {
         for (k = 0; k < NCOLOR; k++) {
            for (l = 0; l < NCOLOR; l++) {
               sum=0.0+0.0*I;
               for (m = 0; m < NCOLOR; m++) {
                  for (n = 0; n < NCOLOR; n++) {
                     sum += B->comp[i][n][m][l] * C->comp[n][j][k][m];
                  }
               }
               A->comp[i][j][k][l] = sum;
            }
         }
      }
   }
}

inline double retr_TensProd(TensProd const * const restrict A)
  {
  int i0, i1;
  double complex tr;
  double ris;

  tr=0.0;
  for(i0=0; i0<NCOLOR; i0++)
     {
     for(i1=0; i1<NCOLOR; i1++)
        {
        tr+=A->comp[i0][i0][i1][i1];
        }
     }
  ris=creal(tr);
  ris/=(NCOLOR*NCOLOR);

  return ris;
  }


inline double imtr_TensProd(TensProd const * const restrict A)
  {
  int i0, i1;
  double complex tr;
  double ris;

  tr=0.0;
  for(i0=0; i0<NCOLOR; i0++)
     {
     for(i1=0; i1<NCOLOR; i1++)
        {
        tr+=A->comp[i0][i0][i1][i1];
        }
     }
  ris=cimag(tr);
  ris/=(NCOLOR*NCOLOR);

  return ris;
  }

void LU_TensProd(TensProd const * const restrict TP, TensProd * restrict result, int * restrict sign);
inline complex double det_TensProd(TensProd const * const TP) {
#ifdef __INTEL_COMPILER
  __assume_aligned(&(TP->comp), DOUBLE_ALIGN);
#endif
   int sign = 0;
   TensProd lu; LU_TensProd(TP, &lu, &sign);

   complex double result = sign > 0 ? 1. + I*0. : -1. + I*0.;
   for (int i = 0; i < NCOLOR; i++) for (int j = 0; j < NCOLOR; j++) {
      result *= lu.comp[i][i][j][j]; // [in1][out1][out2][in2], diagonal is (in1, in2) == (out1, out2)
   }
   return result;
}

void print_on_screen_TensProd(TensProd const * const A);
void print_on_file_TensProd(FILE *fp, TensProd const * const A);
void print_on_binary_file_noswap_TensProd(FILE *fp, TensProd const * const A);
void print_on_binary_file_swap_TensProd(FILE *fp, TensProd const * const A);
void print_on_binary_file_bigen_TensProd(FILE *fp, TensProd const * const A);
void read_from_file_TensProd(FILE *fp, TensProd *A);
void read_from_binary_file_noswap_TensProd(FILE *fp, TensProd *A);
void read_from_binary_file_swap_TensProd(FILE *fp, TensProd *A);
void read_from_binary_file_bigen_TensProd(FILE *fp, TensProd *A);

#endif // TENS_PROD_H

