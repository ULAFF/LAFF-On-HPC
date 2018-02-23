#include <stdio.h>
#include <stdlib.h>

#define alpha( i,j ) A[ (j)*ldA + (i) ]   // map alpha( i,j ) to array A
#define beta( i,j )  B[ (j)*ldB + (i) ]   // map beta( i,j ) to array B
#define gamma( i,j ) C[ (j)*ldC + (i) ]   // map gamma( i,j ) to array C

#define MR 4
#define NR 4
#define KR 4

void Gemm_IJP_MRxNRxKR( int, int, int, double *, int, double *, int, double *, int );

void Gemm_MRxNRxKR( double *, int, double *, int, double *, int );

void GemmWrapper( int m, int n, int k, double *A, int ldA,
		  double *B, int ldB, double *C, int ldC )
{
  if ( m % MR != 0 || n % NR != 0 || k % KR != 0 ){
    printf( "m, n, and k must be multiples of MR, NR, and KR, respectively \n" );
    exit( 0 );
  }
  
  Gemm_IJP_MRxNRxKR( m, n, k, A, ldA, B, ldB, C, ldC );
}

void Gemm_IJP_MRxNRxKR( int m, int n, int k, double *A, int ldA,
                           double *B, int ldB, double *C, int ldC )
{
  for ( int i=0; i<m; i+=MR ) /* m is assumed to be a multiple of MR */
    for ( int j=0; j<n; j+=NR ) /* n is assumed to be a multiple of NR */
      for ( int p=0; p<k; p+=KR ) /* k is assumed to be a multiple of KR */
        Gemm_MRxNRxKR( &alpha( i,p ), ldA, &beta( p,j ), ldB, &gamma( i,j ), ldC );
}

void Gemm_MRxNRxKR( double *A, int ldA, double *B, int ldB, double *C, int ldC )
{
  for ( int i=0; i<MR; i++ )
    for ( int j=0; j<NR; j++ )
      for ( int p=0; p<KR; p++ )
        gamma( i,j ) += alpha( i,p ) * beta( p,j );
}
