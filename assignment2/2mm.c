/**
 * 2mm.c: This file is part of the PolyBench/C 3.2 test suite.
 *
 *
 * Contact: Louis-Noel Pouchet <pouchet@cse.ohio-state.edu>
 * Web address: http://polybench.sourceforge.net
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#ifdef CONFIG_SSE2
#include <emmintrin.h>
#endif

  static
void init_array(int ni, int nj, int nk, int nl,
    double *alpha,
    double *beta,
    double A[ni][nl],
    double B[nk][nj],
    double C[nl][nj],
    double D[ni][nl])
{
  int i, j;

  *alpha = 32412;
  *beta = 2123;
  for (i = 0; i < ni; i++)
    for (j = 0; j < nk; j++)
      A[i][j] = ((double) i*j) / ni;
#ifdef CONFIG_VECTOR
 for (i = 0; i < nk; i++)
    for (j = 0; j < nj; j++)
      B[j][i] = ((double) i*(j+1)) / nj;
#else
  for (i = 0; i < nk; i++)
    for (j = 0; j < nj; j++)
      B[i][j] = ((double) i*(j+1)) / nj;
#endif

#ifdef CONFIG_VECTOR
  for (i = 0; i < nl; i++)
    for (j = 0; j < nj; j++)
      C[j][i] = ((double) i*(j+3)) / nl;
#else
  for (i = 0; i < nl; i++)
    for (j = 0; j < nj; j++)
      C[i][j] = ((double) i*(j+3)) / nl;
#endif

  for (i = 0; i < ni; i++)
    for (j = 0; j < nl; j++)
      D[i][j] = ((double) i*(j+2)) / nk;
}




  static
void print_array(int ni, int nl,
    double D[ni][nl])
{
  int i, j;

  for (i = 0; i < ni; i++)
  {
    for (j = 0; j < nl; j++) {
      fprintf (stderr, "%0.2lf ", D[i][j]);
/*      if ((i * ni + j) % 20 == 0) fprintf (stderr, "\n");*/
    }
    fprintf(stderr,"\n");
    
  }
  fprintf (stderr, "\n");
}




  static
void kernel_2mm(int ni, int nj, int nk, int nl,
    double alpha,
    double beta,
    double tmp[ni][nj],
    double A[ni][nk],
    double B[nk][nj],
    double C[nl][nj],
    double D[ni][nl])
{
  int i, j, k;

#pragma scop

  for (i = 0; i < ni; i++)
    for (j = 0; j < nj; j++)
    {
      tmp[i][j] = 0;
#ifdef CONFIG_SSE2
      __m128d v_subsum = _mm_set1_pd(0);
      __m128d v_alpha = _mm_set1_pd(alpha);
      double subsum[2];
      for (k = 0; k < nk; k +=2)
      {
        //tmp[i][j] += alpha * A[i][k] * B[k][j];
        __m128d v_a = _mm_load_pd(&A[i][k]);
        __m128d v_b = _mm_load_pd(&B[j][k]);
        __m128d v_tmp;
        v_tmp = _mm_mul_pd(v_a, v_b);
        v_tmp = _mm_mul_pd(v_tmp,v_alpha);
        v_subsum = _mm_add_pd(v_subsum,v_tmp);
      }
      _mm_store_pd(subsum,v_subsum);
      tmp[i][j] = subsum[0] + subsum[1];
#else
      for (k = 0; k < nk; ++k)
        tmp[i][j] += alpha * A[i][k] * B[k][j];
#endif
    }


  for (i = 0; i < ni; i++)
    for (j = 0; j < nl; j++)
    {      
      D[i][j] *= beta;
#ifdef CONFIG_SSE2
      __m128d v_subsum = _mm_set1_pd(0);
      double subsum[2];
      for (k = 0; k < nj; k += 2)
      {
        __m128d v_a = _mm_load_pd(&tmp[i][k]);
        __m128d v_c = _mm_load_pd(&C[j][k]);
        __m128d v_tmp;
        v_tmp = _mm_mul_pd(v_a, v_c);
        v_subsum = _mm_add_pd(v_subsum,v_tmp);
        //D[i][j] += tmp[i][k] * C[k][j];
      }
      _mm_store_pd(subsum,v_subsum);
      D[i][j] +=  (subsum[0] + subsum[1]);
#else
      for (k = 0; k < nj; ++k)
        D[i][j] += tmp[i][k] * C[k][j];
#endif      
    }
#pragma endscop

}


int main(int argc, char** argv)
{
  int dump_code = atoi(argv[1]);
  int ni = atoi(argv[2]);
  int nj = atoi(argv[3]);
  int nk = atoi(argv[4]);
  int nl = atoi(argv[5]);


  double alpha;
  double beta;
  double (*tmp)[ni][nj]; tmp = (double(*)[ni][nj])malloc((ni) * (nj) * sizeof(double));;
  double (*A)[ni][nk]; A = (double(*)[ni][nk])malloc((ni) * (nk) * sizeof(double));;
  double (*B)[nk][nj]; B = (double(*)[nk][nj])malloc((nk) * (nj) * sizeof(double));;
  double (*C)[nl][nj]; C = (double(*)[nl][nj])malloc((nl) * (nj) * sizeof(double));;
  double (*D)[ni][nl]; D = (double(*)[ni][nl])malloc((ni) * (nl) * sizeof(double));;


  init_array (ni, nj, nk, nl, &alpha, &beta,
      *A,
      *B,
      *C,
      *D);



  kernel_2mm (ni, nj, nk, nl,
      alpha, beta,
      *tmp,
      *A,
      *B,
      *C,
      *D);





  if (dump_code == 1) print_array(ni, nl, *D);


  free((void*)tmp);;
  free((void*)A);;
  free((void*)B);;
  free((void*)C);;
  free((void*)D);;

  return 0;
}
