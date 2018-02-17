#define alpha( i,j ) A[ (j)*ldA + i ]   // map alpha( i,j ) to array A 
#define beta( i,j )  B[ (j)*ldB + i ]   // map beta( i,j )  to array B
#define gamma( i,j ) C[ (j)*ldC + i ]   // map gamma( i,j ) to array C

void Gemm_J_Gemv_J_Axpy( int, int, int, double *, int, double *, int, double *, int );

void Gemv_J_Axpy( int, int, double *, int, double *, int, double *, int );

void GemmWrapper( int m, int n, int k, double *A, int ldA,
		  double *B, int ldB, double *C, int ldC )
{
  Gemm_J_Gemv_J_Axpy( m, n, k, A, ldA, B, ldB, C, ldC );
}


void Gemm_J_Gemv_J_Axpy( int m, int n, int k,
		   double *A, int ldA,
		   double *B, int ldB,
		   double *C, int ldC )
{
  for ( int j=0; j<n; j++ )
    Gemv_J_Axpy(  ,  ,  ,  ,  ,  ,  ,  );
}
  
