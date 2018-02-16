#define alpha( i,j ) A[ (j)*ldA + i ]   // map alpha( i,j ) to array A 
#define beta( i,j )  B[ (j)*ldB + i ]   // map beta( i,j )  to array B
#define gamma( i,j ) C[ (j)*ldC + i ]   // map gamma( i,j ) to array C

void Gemm_J_Gemv_I_Dots( int, int, int, double *, int, double *, int, double *, int );


void Gemv_I_Dots( int, int, double *, int, double *, int, double *, int );

void GemmWrapper( int m, int n, int k, double *A, int ldA,
		  double *B, int ldB, double *C, int ldC )
{
  Gemm_J_Gemv_I_Dots( m, n, k, A, ldA, B, ldB, C, ldC );
}


void Gemm_J_Gemv_I_Dots( int m, int n, int k,
		   double *A, int ldA,
		   double *B, int ldB,
		   double *C, int ldC )
{
  for ( int j=0; j<n; j++ )
    Gemv_I_Dots( m, k, A, ldA, &beta( 0, j ), 1, &gamma( 0,j ), 1 );
}
  
