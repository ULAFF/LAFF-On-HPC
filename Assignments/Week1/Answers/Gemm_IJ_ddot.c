#define alpha( i,j ) A[ (j)*ldA + i ]   // map alpha( i,j ) to array A 
#define beta( i,j )  B[ (j)*ldB + i ]   // map beta( i,j )  to array B
#define gamma( i,j ) C[ (j)*ldC + i ]   // map gamma( i,j ) to array C

void Gemm_IJ_ddot( int, int, int, const double *, int, const double *, int, double *, int );
double ddot_( int *, const double *, int *, const double *, int * );

void GemmWrapper( int m, int n, int k, const double *A, int ldA,
		  const double *B, int ldB, double *C, int ldC )
{
  Gemm_IJ_ddot( m, n, k, A, ldA, B, ldB, C, ldC );
}


void Gemm_IJ_ddot( int m, int n, int k, const double *A, int ldA,
                   const double *B, int ldB, double *C, int ldC )
{
  for ( int i=0; i<m; i++ )
    for ( int j=0; j<n; j++ ){
	int i_one = 1;      // declared to be able to pass "1" by address 
	gamma( i,j ) +=  ddot_( &k, &alpha( i,0 ), &ldA, 
                                    &beta( 0,j ), &i_one );
      }
}
