#define alpha( i,j ) A[ (j)*ldA + i ]   // map alpha( i,j ) to array A 
#define chi( i )  x[ (i)*incx ]         // map chi( i )  to array x
#define psi( i )  y[ (i)*incy ]         // map psi( i )  to array x

void Dots( int, double *, int, double *, int, double * );

void Gemv_I_Dots( int m, int n, double *A, int ldA,
                  double *x, int incx, double *y, int incy )
{
  for ( int i=0; i<m; i++ )
    Dots( n, &alpha( i,0 ), ldA, x, incx, &psi( i ) );
}
