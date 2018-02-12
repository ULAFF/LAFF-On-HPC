#include <blis.h>

#define alpha( i,j ) A[ (j)*ldA + i ]   // map alpha( i,j ) to array A 
#define beta( i,j )  B[ (j)*ldB + i ]   // map beta( i,j )  to array B
#define gamma( i,j ) C[ (j)*ldC + i ]   // map gamma( i,j ) to array C

void Gemm_JP_bli_daxpyv( int, int, int, double *, int, double *, int, double *, int );

void GemmWrapper( int m, int n, int k, double *A, int ldA,
		  double *B, int ldB, double *C, int ldC )
{
  Gemm_JP_bli_daxpyv( m, n, k, A, ldA, B, ldB, C, ldC );
}


void Gemm_JP_bli_daxpyv( int m, int n, int k,
			 double *A, int ldA,
			 double *B, int ldB,
			 double *C, int ldC )
{
  for ( int j=0; j<n; j++ )
    for ( int p=0; p<k; p++ )
      bli_daxpyv( BLIS_NO_CONJUGATE, n, &beta( p,j ), &alpha( 0,p ), 1, &gamma( 0,j ), 1, NULL );
}
  
