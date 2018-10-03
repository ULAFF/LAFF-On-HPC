#define alpha( i,j ) A[ (j)*ldA + i ]   // map alpha( i,j ) to array A 
#define chi( i )  x[ (i)*incx ]         // map chi( i )  to array x
#define psi( i )  y[ (i)*incy ]         // map psi( i )  to array x

void Axpy( int, double, double *, int, double *, int );

void Ger( int m, int n, double *x, int incx,
		 double *y, int incy, double *A, int ldA )
{
  for ( int i=0; i<m; i++ )
    Axpy( n, chi( i ), y, incy, &alpha( i,0 ), ldA );
}
