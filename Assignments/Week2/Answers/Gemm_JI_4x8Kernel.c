#include <stdio.h>
#include <stdlib.h>

#define alpha( i,j ) A[ (j)*ldA + (i) ]   // map alpha( i,j ) to array A
#define beta( i,j )  B[ (j)*ldB + (i) ]   // map beta( i,j ) to array B
#define gamma( i,j ) C[ (j)*ldC + (i) ]   // map gamma( i,j ) to array C

#define MR 4
#define NR 8

void Gemm_JI_4x8Kernel( int, int, int, double *, int, double *, int, double *, int );

void Gemm_4x8Kernel( int, double *, int, double *, int, double *, int );

void GemmWrapper( int m, int n, int k, double *A, int ldA,
		  double *B, int ldB, double *C, int ldC )
{
  if ( m % MR != 0 || n % NR != 0 ){
    printf( "m and n must be multiples of MR and NR, respectively \n" );
    exit( 0 );
  }
  
  Gemm_JI_4x8Kernel( m, n, k, A, ldA, B, ldB, C, ldC );
}

void Gemm_JI_4x8Kernel( int m, int n, int k, double *A, int ldA,
                                  double *B, int ldB, double *C, int ldC )
{
  for ( int j=0; j<n; j+=NR ) /* n is assumed to be a multiple of NR */
    for ( int i=0; i<m; i+=MR ) /* m is assumed to be a multiple of MR */
      Gemm_4x8Kernel( k, &alpha( i,0 ), ldA, &beta( 0,j ), ldB, &gamma( i,j ), ldC );
}

#include<immintrin.h>

void Gemm_4x8Kernel( int k, double *A, int ldA, double *B, int ldB,
		double *C, int ldC )
{
  /* Declare vector registers to hold 4x8 C and load them */
  __m256d gamma_0123_0 = _mm256_loadu_pd( &gamma( 0,0 ) );
  __m256d gamma_0123_1 = _mm256_loadu_pd( &gamma( 0,1 ) );
  __m256d gamma_0123_2 = _mm256_loadu_pd( &gamma( 0,2 ) );
  __m256d gamma_0123_3 = _mm256_loadu_pd( &gamma( 0,3 ) );
  __m256d gamma_0123_4 = _mm256_loadu_pd( &gamma( 0,4 ) );
  __m256d gamma_0123_5 = _mm256_loadu_pd( &gamma( 0,5 ) );
  __m256d gamma_0123_6 = _mm256_loadu_pd( &gamma( 0,6 ) );
  __m256d gamma_0123_7 = _mm256_loadu_pd( &gamma( 0,7 ) );
   	
  for ( int p=0; p<k; p++ ){
    /* Declare vector register for load/broadcasting beta( b,j ) */
    __m256d beta_p_j;
    
    /* Declare a vector register to hold the current column of A and load
       it with the four elements of that column. */
    __m256d alpha_0123_p = _mm256_loadu_pd( &alpha( 0,p ) );

    /* Load/broadcast beta( p,0 ). */
    beta_p_j = _mm256_broadcast_sd( &beta( p, 0) );
    
    /* update the first column of C with the current column of A times
       beta ( p,0 ) */
    gamma_0123_0 = _mm256_fmadd_pd( alpha_0123_p, beta_p_j, gamma_0123_0 );
    
    /* REPEAT for second, third, and fourth columns of C.  Notice that the 
       current column of A needs not be reloaded. */

    /* Load/broadcast beta( p,1 ). */
    beta_p_j = _mm256_broadcast_sd( &beta( p, 1) );
    
    /* update the second column of C with the current column of A times
       beta ( p,1 ) */
    gamma_0123_1 = _mm256_fmadd_pd( alpha_0123_p, beta_p_j, gamma_0123_1 );

    /* Load/broadcast beta( p,2 ). */
    beta_p_j = _mm256_broadcast_sd( &beta( p, 2) );
    
    /* update the third column of C with the current column of A times
       beta ( p,2 ) */
    gamma_0123_2 = _mm256_fmadd_pd( alpha_0123_p, beta_p_j, gamma_0123_2 );

    /* Load/broadcast beta( p,3 ). */
    beta_p_j = _mm256_broadcast_sd( &beta( p, 3) );
    
    /* update the fourth column of C with the current column of A times
       beta ( p,3 ) */
    gamma_0123_3 = _mm256_fmadd_pd( alpha_0123_p, beta_p_j, gamma_0123_3 );

    /* Load/broadcast beta( p,4 ). */
    beta_p_j = _mm256_broadcast_sd( &beta( p, 4) );
    
    /* update the fifth column of C with the current column of A times
       beta ( p,4 ) */
    gamma_0123_4 = _mm256_fmadd_pd( alpha_0123_p, beta_p_j, gamma_0123_4 );

    /* Load/broadcast beta( p,5 ). */
    beta_p_j = _mm256_broadcast_sd( &beta( p, 5) );
    
    /* update the sixth column of C with the current column of A times
       beta ( p,5 ) */
    gamma_0123_5 = _mm256_fmadd_pd( alpha_0123_p, beta_p_j, gamma_0123_5 );

    /* Load/broadcast beta( p,6 ). */
    beta_p_j = _mm256_broadcast_sd( &beta( p, 6) );
    
    /* update the seventh column of C with the current column of A times
       beta ( p,6 ) */
    gamma_0123_6 = _mm256_fmadd_pd( alpha_0123_p, beta_p_j, gamma_0123_6 );
    
    /* Load/broadcast beta( p,7 ). */
    beta_p_j = _mm256_broadcast_sd( &beta( p, 7) );
    
    /* update the eighth column of C with the current column of A times
       beta ( p,7 ) */
    gamma_0123_7 = _mm256_fmadd_pd( alpha_0123_p, beta_p_j, gamma_0123_7 );
   
  }
  
  /* Store the updated results */
  _mm256_storeu_pd( &gamma(0,0), gamma_0123_0 );
  _mm256_storeu_pd( &gamma(0,1), gamma_0123_1 );
  _mm256_storeu_pd( &gamma(0,2), gamma_0123_2 );
  _mm256_storeu_pd( &gamma(0,3), gamma_0123_3 );
  _mm256_storeu_pd( &gamma(0,4), gamma_0123_4 );
  _mm256_storeu_pd( &gamma(0,5), gamma_0123_5 );
  _mm256_storeu_pd( &gamma(0,6), gamma_0123_6 );
  _mm256_storeu_pd( &gamma(0,7), gamma_0123_7 );
}
