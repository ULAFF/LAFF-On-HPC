#include <stdio.h>
#include <stdlib.h>

#define alpha( i,j ) A[ (j)*ldA + (i) ]   // map alpha( i,j ) to array A
#define beta( i,j )  B[ (j)*ldB + (i) ]   // map beta( i,j ) to array B
#define gamma( i,j ) C[ (j)*ldC + (i) ]   // map gamma( i,j ) to array C

#define min( x, y ) ( (x) < (y) ? (x) : (y) )

#define MB 4
#define NB 4
#define KB 4

void Gemm_mbxnbxkb( int, int, int, double *, int, double *, int, double *, int );

void MyGemm( int m, int n, int k, double *A, int ldA,
	     double *B, int ldB, double *C, int ldC )
{
  for ( int i=0; i<m; i+=MB ){
    int ib = min( m-i, MB );    /* Size for "finge" block */ 
    for ( int j=0; j<n; j+=NB ){
      int jb = min( n-j, NB );    /* Size for "finge" block */ 
      for ( int p=0; p<k; p+=KB ){ 
        int pb = min( k-p, KB );    /* Size for "finge" block */ 
        Gemm_mbxnbxkb( ib, jb, pb, &alpha( i,p ), ldA, &beta( p,j ), ldB,
		                   &gamma( i,j ), ldC );
      }
    }
  }
}

void Gemm_mbxnbxkb( int m, int n, int k, double *A, int ldA, 
		    double *B, int ldB,  double *C, int ldC )
{
  for ( int p=0; p<k; p++ )
    for ( int j=0; j<n; j++ )
      for ( int i=0; i<m; i++ )
        gamma( i,j ) += alpha( i,p ) * beta( p,j );
}
