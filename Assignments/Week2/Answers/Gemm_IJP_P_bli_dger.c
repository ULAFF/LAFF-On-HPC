#include <stdio.h>
#include <blis.h>

#define alpha( i,j ) A[ (j)*ldA + (i) ]   // map alpha( i,j ) to array A
#define beta( i,j )  B[ (j)*ldB + (i) ]   // map beta( i,j ) to array B
#define gamma( i,j ) C[ (j)*ldC + (i) ]   // map gamma( i,j ) to array C

#define MR 4
#define NR 4
#define KR 4

void Gemm_IJP_P_bli_dger( int, int, int, double *, int, double *, int, double *, int );

void Gemm_P_bli_dger( int, int, int, double *, int, double *, int, double *, int );

void GemmWrapper( int m, int n, int k, double *A, int ldA,
		  double *B, int ldB, double *C, int ldC )
{
  if ( m % MR != 0 || n % NR != 0 || k % KR != 0 ){
    printf( "m, n, and k must be multiples of MR, NR, and KR, respectively \n" );
    exit( 0 );
  }
  
  Gemm_IJP_P_bli_dger( m, n, k, A, ldA, B, ldB, C, ldC );
}

void Gemm_IJP_P_bli_dger( int m, int n, int k, double *A, int ldA,
			  double *B, int ldB, double *C, int ldC )
{
  for ( int i=0; i<m; i+=MR ) /* m is assumed to be a multiple of MR */
    for ( int j=0; j<n; j+=NR ) /* n is assumed to be a multiple of NR */
      for ( int p=0; p<k; p+=KR ) /* k is assumed to be a multiple of KR */
	Gemm_P_bli_dger( MR, NR, KR, &alpha( i,p ), ldA,
			 &beta( p,j ), ldB, &gamma( i,j ), ldC );
}

void Gemm_P_bli_dger( int m, int n, int k, double *A, int ldA,
		      double *B, int ldB, double *C, int ldC )
{
  double d_one = 1.0;
  
  for ( int p=0; p<k; p++ )
    bli_dger( BLIS_NO_CONJUGATE, BLIS_NO_CONJUGATE, m, n, &d_one, &alpha( 0,p ), 1,
              &beta( p,0 ), ldB, C, 1, ldC, NULL );
}
  
