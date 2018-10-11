#include <stdio.h>
#include <stdlib.h>

#define alpha( i,j ) A[ (j)*ldA + (i) ]   // map alpha( i,j ) to array A
#define beta( i,j )  B[ (j)*ldB + (i) ]   // map beta( i,j ) to array B
#define gamma( i,j ) C[ (j)*ldC + (i) ]   // map gamma( i,j ) to array C

#define min( x, y ) ( ( x ) < ( y ) ? x : y )

#define NC 1024
#define KC 128
#define MC 180
#define MR 4
#define NR 4

void LoopFive( int, int, int, double *, int, double *, int, double *, int );
void LoopFour( int, int, int, double *, int, double *, int, double *, int );
void LoopThree( int, int, int, double *, int, double *, int, double *, int );
void LoopTwo( int, int, int, double *, int, double *, int, double *, int );
void LoopOne( int, int, int, double *, int, double *, int, double *, int );
void Gemm_4x4Kernel( int, double *, int, double *, int, double *, int );

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
    Gemm_4x4Kernel( k, &alpha( i, 0 ), ldA, B, ldB, &gamma( i,0 ), ldC );
  }
}

#include<immintrin.h>

void Gemm_4x4Kernel( int k, double *A, int ldA, double *B, int ldB,
		double *C, int ldC )
{
  /* Declare vector registers to hold 4x4 C and load them */
  __m256d gamma_0123_0 = _mm256_loadu_pd( &gamma( 0,0 ) );
  __m256d gamma_0123_1 = _mm256_loadu_pd( &gamma( 0,1 ) );
  __m256d gamma_0123_2 = _mm256_loadu_pd( &gamma( 0,2 ) );
  __m256d gamma_0123_3 = _mm256_loadu_pd( &gamma( 0,3 ) );
   	
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
  }
  
  /* Store the updated results */
  _mm256_storeu_pd( &gamma(0,0), gamma_0123_0 );
  _mm256_storeu_pd( &gamma(0,1), gamma_0123_1 );
  _mm256_storeu_pd( &gamma(0,2), gamma_0123_2 );
  _mm256_storeu_pd( &gamma(0,3), gamma_0123_3 );
}
