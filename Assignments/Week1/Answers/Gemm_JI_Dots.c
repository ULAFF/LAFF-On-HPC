#define alpha( i,j ) A[ (j)*ldA + i ]   // map alpha( i,j ) to array A 
#define beta( i,j )  B[ (j)*ldB + i ]   // map beta( i,j )  to array B
#define gamma( i,j ) C[ (j)*ldC + i ]   // map gamma( i,j ) to array C

void Gemm_JI_Dots( int, int, int, double *, int, double *, int, double *, int );
void Dots( int, double *, int, double *, int, double * );

void GemmWrapper( int m, int n, int k, double *A, int ldA,
		  double *B, int ldB, double *C, int ldC )
{
  Gemm_JI_Dots( m, n, k, A, ldA, B, ldB, C, ldC );
}


void Gemm_JI_Dots( int m, int n, int k, double *A, int ldA,
                   double *B, int ldB, double *C, int ldC )
{
  for ( int j=0; j<n; j++ )
    for ( int i=0; i<m; i++ )
      Dots( k, &alpha( i,0 ), ldA, &beta( 0,j ), 1, &gamma( i,j ) );
}
