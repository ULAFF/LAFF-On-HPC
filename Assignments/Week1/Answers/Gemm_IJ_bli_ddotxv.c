#include "blis.h"

#define alpha( i,j ) A[ (j)*ldA + i ]   // map alpha( i,j ) to array A 
#define beta( i,j )  B[ (j)*ldB + i ]   // map beta( i,j )  to array B
#define gamma( i,j ) C[ (j)*ldC + i ]   // map gamma( i,j ) to array C

void Gemm_IJ_bli_ddotxv( int, int, int, const double *, int, const double *, int, double *, int );

void GemmWrapper( int m, int n, int k, const double *A, int ldA,
		  const double *B, int ldB, double *C, int ldC )
{
  Gemm_IJ_bli_ddotxv( m, n, k, A, ldA, B, ldB, C, ldC );
}


void Gemm_IJ_bli_ddotxv( int m, int n, int k, const double *A, int ldA,
                   const double *B, int ldB, double *C, int ldC )
{
  for ( int i=0; i<m; i++ )
    for ( int j=0; j<n; j++ ){
      double d_one=1.0;   // Needed to pass 1.0 by address
      bli_ddotxv( BLIS_NO_CONJUGATE, BLIS_NO_CONJUGATE,
                  k, &d_one, &alpha( i,0 ), ldA, &beta( 0,j ), 1,
                  &d_one, &gamma( i,j ), NULL );
    }
}
