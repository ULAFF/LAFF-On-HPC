#define chi( i ) x[ (i)*incx ]   // map chi( i ) to array x 
#define psi( i ) y[ (i)*incy ]   // map psi( i ) to array y

void Dots( int n, const double *x, int incx, const double *y, int incy, double *gamma )
{
  for ( int i=0; i<n; i++ )
    *gamma += chi( i ) * psi( i );
}
