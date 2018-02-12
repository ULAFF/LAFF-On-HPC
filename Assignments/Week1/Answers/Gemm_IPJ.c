#define alpha( i,j ) A[ (j)*ldA + i ]   // map alpha( i,j ) to array A 
#define beta( i,j )  B[ (j)*ldB + i ]   // map beta( i,j )  to array B
#define gamma( i,j ) C[ (j)*ldC + i ]   // map gamma( i,j ) to array C

void Gemm_IPJ( int, int, int, const double *, int, const double *, int, double *, int );

void GemmWrapper( int m, int n, int k, const double *A, int ldA,
		  const double *B, int ldB, double *C, int ldC )
{
  Gemm_IPJ( m, n, k, A, ldA, B, ldB, C, ldC );
}


void Gemm_IPJ( int m, int n, int k,
               const double *A, int ldA,
               const double *B, int ldB,
               double *C, int ldC )
{
  for ( int i=0; i<m; i++ )
    for ( int p=0; p<k; p++ )
      for ( int j=0; j<n; j++ )
        gamma( i,j ) += alpha( i,p ) * beta( p,j );
}
  
