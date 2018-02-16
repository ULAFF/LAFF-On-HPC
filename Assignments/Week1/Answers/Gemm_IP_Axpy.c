#define alpha( i,j ) A[ (j)*ldA + i ]   // map alpha( i,j ) to array A 
#define beta( i,j )  B[ (j)*ldB + i ]   // map beta( i,j )  to array B
#define gamma( i,j ) C[ (j)*ldC + i ]   // map gamma( i,j ) to array C

void Gemm_IP_Axpy( int, int, int, double *, int, double *, int, double *, int );
void Axpy( int, double, double *, int, double *, int );

void GemmWrapper( int m, int n, int k, double *A, int ldA,
		  double *B, int ldB, double *C, int ldC )
{
  Gemm_IP_Axpy( m, n, k, A, ldA, B, ldB, C, ldC );
}


void Gemm_IP_Axpy( int m, int n, int k,
		   double *A, int ldA,
		   double *B, int ldB,
		   double *C, int ldC )
{
  for ( int i=0; i<m; i++ )
    for ( int p=0; p<k; p++ )
      Axpy( n, alpha( i,p ), &beta( p,0 ), ldB, &gamma( i,0 ), ldC );
}
  
