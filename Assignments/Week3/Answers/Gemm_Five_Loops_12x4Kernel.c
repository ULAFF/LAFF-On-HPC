#include <stdio.h>
#include <stdlib.h>

#define alpha( i,j ) A[ (j)*ldA + (i) ]   // map alpha( i,j ) to array A
#define beta( i,j )  B[ (j)*ldB + (i) ]   // map beta( i,j ) to array B
#define gamma( i,j ) C[ (j)*ldC + (i) ]   // map gamma( i,j ) to array C

#define min( x, y ) ( ( x ) < ( y ) ? x : y )

#define NC 1024
#define KC 128
#define MC 180
#define MR 12
#define NR 4

void LoopFive( int, int, int, double *, int, double *, int, double *, int );
void LoopFour( int, int, int, double *, int, double *, int, double *, int );
void LoopThree( int, int, int, double *, int, double *, int, double *, int );
void LoopTwo( int, int, int, double *, int, double *, int, double *, int );
void LoopOne( int, int, int, double *, int, double *, int, double *, int );
void Gemm_12x4Kernel( int, double *, int, double *, int, double *, int );
  
void MyGemm( int m, int n, int k, double *A, int ldA,
	     double *B, int ldB, double *C, int ldC )
{
  if ( m % MR != 0 || MC % MR != 0 ){
    printf( "m and MC must be multiples of MR\n" );
    exit( 0 );
  }
  if ( n % NR != 0 || NC % NR != 0 ){
    printf( "n and NC must be multiples of NR\n" );
    exit( 0 );
  }

  LoopFive( m, n, k, A, ldA, B, ldB, C, ldC );
}

void LoopFive( int m, int n, int k, double *A, int ldA,
		   double *B, int ldB, double *C, int ldC )
{
  for ( int j=0; j<n; j+=NC ) {
    int jb = min( NC, n-j );    /* Last loop may not involve a full block */
    LoopFour( m, jb, k, A, ldA, &beta( 0,j ), ldB, &gamma( 0,j ), ldC );
  }
}

void LoopFour( int m, int n, int k, double *A, int ldA, 
	       double *B, int ldB, double *C, int ldC )
{
  for ( int p=0; p<k; p+=KC ) {
    int pb = min( KC, k-p );    /* Last loop may not involve a full block */
    LoopThree( m, n, pb, &alpha( 0, p ), ldA, &beta( p, 0 ), ldB, C, ldC );
  }
}

void LoopThree( int m, int n, int k, double *A, int ldA, 
		double *B, int ldB, double *C, int ldC )
{
  for ( int i=0; i<m; i+=MC ) {
    int ib = min( MC, m-i );    /* Last loop may not involve a full block */
    LoopTwo( ib, n, k, &alpha( i, 0), ldA, B, ldB, &gamma( i,0 ), ldC );
  }
}

void LoopTwo( int m, int n, int k, double *A, int ldA,
	      double *B, int ldB, double *C, int ldC )
{
  for ( int j=0; j<n; j+=NR ) {
    int jb = min( NR, n-j );
    LoopOne( m, jb, k, A, ldA, &beta( 0,j ), ldB, &gamma( 0,j ), ldC );
  }
}

void LoopOne( int m, int n, int k, double *A, int ldA,
	      double *B, int ldB, double *C, int ldC )
{
  for ( int i=0; i<m; i+=MR ) {
    int ib = min( MR, m-i );
    Gemm_12x4Kernel( k, &alpha( i, 0 ), ldA, B, ldB, &gamma( i,0 ), ldC );
  }
}

#include<immintrin.h>

void Gemm_12x4Kernel( int k, double *A, int ldA,
                             double *B, int ldB, double *C, int ldC )
{
  __m256d gamma_0123_0 = _mm256_loadu_pd( &gamma( 0,0 ) );
  __m256d gamma_4567_0 = _mm256_loadu_pd( &gamma( 4,0 ) );
  __m256d gamma_891011_0 = _mm256_loadu_pd( &gamma( 8,0 ) );
  __m256d gamma_0123_1 = _mm256_loadu_pd( &gamma( 0,1 ) );
  __m256d gamma_4567_1 = _mm256_loadu_pd( &gamma( 4,1 ) );
  __m256d gamma_891011_1 = _mm256_loadu_pd( &gamma( 8,1 ) );
  __m256d gamma_0123_2 = _mm256_loadu_pd( &gamma( 0,2 ) );
  __m256d gamma_4567_2 = _mm256_loadu_pd( &gamma( 4,2 ) );
  __m256d gamma_891011_2 = _mm256_loadu_pd( &gamma( 8,2 ) );
  __m256d gamma_0123_3 = _mm256_loadu_pd( &gamma( 0,3 ) );
  __m256d gamma_4567_3 = _mm256_loadu_pd( &gamma( 4,3 ) );
  __m256d gamma_891011_3 = _mm256_loadu_pd( &gamma( 8,3 ) );
   	
  for ( int p=0; p<k; p++ ){
    /* Declare three vector registers to hold the current column of A and load
       it with the twelve elements of that column. */
    __m256d alpha_0123_p = _mm256_loadu_pd( &alpha( 0,p ) );
    __m256d alpha_4567_p = _mm256_loadu_pd( &alpha( 4,p ) );
    __m256d alpha_891011_p = _mm256_loadu_pd( &alpha( 8,p ) );

    /* Declare a vector register to hold elements beta( p,0 ) and duplicate
       it in the register. */
    __m256d beta_p_j = _mm256_broadcast_sd( &beta( p, 0) );
    
    /* update the first column of C with the current column of A times
       beta ( p,0 ) */
    gamma_0123_0 = _mm256_fmadd_pd( alpha_0123_p, beta_p_j, gamma_0123_0 );
    gamma_4567_0 = _mm256_fmadd_pd( alpha_4567_p, beta_p_j, gamma_4567_0 );
    gamma_891011_0 = _mm256_fmadd_pd( alpha_891011_p, beta_p_j, gamma_891011_0 );
    
    /* REPEAT for second, third, and fourth columns of C.  Notice that the 
       current column of A needs not be reloaded and the register in which
       beta( p, 0 ) was loaded can be reused. */

    /* duplicate  beta( p,1 ) in the register. */
    beta_p_j = _mm256_broadcast_sd( &beta( p, 1) );
    
    /* update the first column of C with the current column of A times
       beta ( p,1 ) */
    gamma_0123_1 = _mm256_fmadd_pd( alpha_0123_p, beta_p_j, gamma_0123_1 );
    gamma_4567_1 = _mm256_fmadd_pd( alpha_4567_p, beta_p_j, gamma_4567_1 );
    gamma_891011_1 = _mm256_fmadd_pd( alpha_891011_p, beta_p_j, gamma_891011_1 );

    /* duplicate  beta( p,2 ) in the register. */
    beta_p_j = _mm256_broadcast_sd( &beta( p, 2) );
    
    /* update the first column of C with the current column of A times
       beta ( p,2 ) */
    gamma_0123_2 = _mm256_fmadd_pd( alpha_0123_p, beta_p_j, gamma_0123_2 );
    gamma_4567_2 = _mm256_fmadd_pd( alpha_4567_p, beta_p_j, gamma_4567_2 );
    gamma_891011_2 = _mm256_fmadd_pd( alpha_891011_p, beta_p_j, gamma_891011_2 );

    /* duplicate  beta( p,3 ) in the register. */
    beta_p_j = _mm256_broadcast_sd( &beta( p, 3) );
    
    /* update the first column of C with the current column of A times
       beta ( p,3 ) */
    gamma_0123_3 = _mm256_fmadd_pd( alpha_0123_p, beta_p_j, gamma_0123_3 );
    gamma_4567_3 = _mm256_fmadd_pd( alpha_4567_p, beta_p_j, gamma_4567_3 );
    gamma_891011_3 = _mm256_fmadd_pd( alpha_891011_p, beta_p_j, gamma_891011_3 );
  }
  
  /* Store the updated results */
  _mm256_storeu_pd( &gamma(0,0), gamma_0123_0 );
  _mm256_storeu_pd( &gamma(4,0), gamma_4567_0 );
  _mm256_storeu_pd( &gamma(8,0), gamma_891011_0 );
  _mm256_storeu_pd( &gamma(0,1), gamma_0123_1 );
  _mm256_storeu_pd( &gamma(4,1), gamma_4567_1 );
  _mm256_storeu_pd( &gamma(8,1), gamma_891011_1 );
  _mm256_storeu_pd( &gamma(0,2), gamma_0123_2 );
  _mm256_storeu_pd( &gamma(4,2), gamma_4567_2 );
  _mm256_storeu_pd( &gamma(8,2), gamma_891011_2 );
  _mm256_storeu_pd( &gamma(0,3), gamma_0123_3 );
  _mm256_storeu_pd( &gamma(4,3), gamma_4567_3 );
  _mm256_storeu_pd( &gamma(8,3), gamma_891011_3 );

  return;
}
