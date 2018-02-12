#define alpha( i,j ) A[ (j)*ldA + i ]   // map alpha( i,j ) to array A 
#define beta( i,j )  B[ (j)*ldB + i ]   // map beta( i,j )  to array B
#define gamma( i,j ) C[ (j)*ldC + i ]   // map gamma( i,j ) to array C

void Gemm_IP_daxpy( int, int, int, const double *, int, const double *, int, double *, int );
void daxpy_( int *, const double *, const double *, int *, double *, int * );

void GemmWrapper( int m, int n, int k, const double *A, int ldA,
		  const double *B, int ldB, double *C, int ldC )
{
  Gemm_IP_daxpy( m, n, k, A, ldA, B, ldB, C, ldC );
}


void Gemm_IP_daxpy( int m, int n, int k,
		   const double *A, int ldA,
		   const double *B, int ldB,
		   double *C, int ldC )
{
  for ( int i=0; i<m; i++ )
    for ( int p=0; p<k; p++ )
      daxpy_( &n, &alpha( i,p ), &beta( p,0 ), &ldB, &gamma( i,0 ), &ldC );
}
  
