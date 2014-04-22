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

#ifdef CONFIG_TILE
#ifndef N
#define N 2


static void init_tile(double *matrix, int block_row_size, int block_col_size)
{
  memset((void*)matrix, 0, sizeof(double) * block_row_size*block_col_size);
}

//C = scalar * A * B
inline void calculate_first_phase_block(double *A, double *B, double *C, double scalar, int rowsize_A, int colsize_A, int colsize_B)
{
  int i,j,k;
  for(i = 0 ; i < rowsize_A ; i++)
  {
    for(j = 0 ; j < colsize_A ; j++)
    {
      for(k = 0 ; k < colsize_A ; k++)
      {
        *(C + i * rowsize_A + j) = *(A + i * rowsize_A + k) + *(B + colsize_A * k + j) * scalar;
      }
    }
  }
}
#endif
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
  for (i = 0; i < nk; i++)
    for (j = 0; j < nj; j++)
      B[i][j] = ((double) i*(j+1)) / nj;
  for (i = 0; i < nl; i++)
    for (j = 0; j < nj; j++)
      C[i][j] = ((double) i*(j+3)) / nl;
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
    for (j = 0; j < nl; j++) {
      fprintf (stderr, "%0.2lf ", D[i][j]);
      if ((i * ni + j) % 20 == 0) fprintf (stderr, "\n");
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

#ifdef CONFIG_TILE
  int block_ni = ni / N;
  int block_nj = nj / N;
  int block_nk = nk / N;
  int block_nl = nl / N;
#endif
ß
#pragma scop

#ifndef CONFIG_TILE
  for (i = 0; i < ni; i++)
    for (j = 0; j < nj; j++)
    {
      tmp[i][j] = 0;
      for (k = 0; k < nk; ++k)
        tmp[i][j] += alpha * A[i][k] * B[k][j];
    }
  for (i = 0; i < ni; i++)
    for (j = 0; j < nl; j++)
    {
      D[i][j] *= beta;
      for (k = 0; k < nj; ++k)
        D[i][j] += tmp[i][k] * C[k][j];
    }
#else  //Tiling

for(i = 0 ; i < N ; i++)
{
  for(j = 0 ; j < N ; j++)
  {
    init_tile(&tmp[i*block_ni][j*block_nj], block_ni, block_nj);
    for(k = 0 ; k < N ; k++)
    {
      //Multiply blocks
/*
3
1 1 1 1  2 2 2 2  24 24 24 24
1 1 1 1  2 2 2 2  24 24 24 24
1 1 1 1  2 2 2 2  24 24 24 24
1 1 1 1  2 2 2 2  24 24 24 24

0 0
0 0
12 12 12 12  24 24
12 12 12 12  24 24
*/
    }
  }
}
/*for (i = 0; i < ni; i++)
  for (j = 0; j < nj; j++)
  {
    tmp[i][j] = 0;
    for (k = 0; k < nk; ++k)
      tmp[i][j] += alpha * A[i][k] * B[k][j];
  }
for (i = 0; i < ni; i++)
  for (j = 0; j < nl; j++)
  {
    D[i][j] *= beta;
    for (k = 0; k < nj; ++k)
      D[i][j] += tmp[i][k] * C[k][j];
  }*/
#endif
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
