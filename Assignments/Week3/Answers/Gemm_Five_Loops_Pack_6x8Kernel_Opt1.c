#include <stdio.h>
#include <stdlib.h>

#define alpha( j,i ) A[ (j)*ldA + (i) ]   // map alpha( i,j ) to array A
#define beta( j,i )  B[ (j)*ldB + (i) ]   // map beta( i,j ) to array B
#define gamma( j,i ) C[ (j)*ldC + (i) ]   // map gamma( i,j ) to array C

#define min( x, y ) ( ( x ) < ( y ) ? x : y )

#define NC 1024
#define KC 128
#define MC 180
#define MR 6
#define NR 8

void LoopFive( int, int, int, double *, int, double *, int, double *, int );
void LoopFour( int, int, int, double *, int, double *, int,  double *, int );
void LoopThree( int, int, int, double *, int, double *, double *, int );
void LoopTwo( int, int, int, double *, double *, double *, int );
void LoopOne( int, int, int, double *, double *, double *, int );
void Gemm_6x8Kernel_Packed( int, double *, double *, double *, int );
void PackBlockA_MCxKC( int, int, double *, int, double * );
void PackPanelB_KCxNC( int, int, double *, int, double * );
  
void MyGemm( int n, int m, int k, double *B, int ldB,
	     double *A, int ldA, double *C, int ldC )
{
  if ( m % MR != 0 || MC % MR != 0 ){
    printf( "m and MC must be multiples of MR\n" );
    exit( 0 );
  }
  if ( n % NR != 0 || NC % NR != 0 ){
    printf( "n and NC must be multiples of NR\n" );
    exit( 0 );
  }

  LoopFive( m, n, k, A, ldA, B, ldB, C, ldC );
}

void LoopFive( int m, int n, int k, double *A, int ldA,
		   double *B, int ldB, double *C, int ldC )
{
  for ( int j=0; j<n; j+=NC ) {
    int jb = min( NC, n-j );    /* Last loop may not involve a full block */
    LoopFour( m, jb, k, A, ldA, &beta( 0,j ), ldB, &gamma( 0,j ), ldC );
  } 
}

void LoopFour( int m, int n, int k, double *A, int ldA, double *B, int ldB,
	       double *C, int ldC )
{
  double *Btilde = ( double * ) malloc( KC * NC * sizeof( double ) );
  
  for ( int p=0; p<k; p+=KC ) {
    int pb = min( KC, k-p );    /* Last loop may not involve a full block */
    PackPanelB_KCxNC( pb, n, &beta( p, 0 ), ldB, Btilde );
    LoopThree( m, n, pb, &alpha( 0, p ), ldA, Btilde, C, ldC );
  }

  free( Btilde); 
}

void LoopThree( int m, int n, int k, double *A, int ldA, double *Btilde, double *C, int ldC )
{
  double *Atilde = ( double * ) malloc( MC * KC * sizeof( double ) );
       
  for ( int i=0; i<m; i+=MC ) {
    int ib = min( MC, m-i );    /* Last loop may not involve a full block */
    PackBlockA_MCxKC( ib, k, &alpha( i, 0 ), ldA, Atilde );
    LoopTwo( ib, n, k, Atilde, Btilde, &gamma( i,0 ), ldC );
  }

  free( Atilde);
}

void LoopTwo( int m, int n, int k, double *Atilde, double *Btilde, double *C, int ldC )
{
  for ( int j=0; j<n; j+=NR ) {
    int jb = min( NR, n-j );
    LoopOne( m, jb, k, Atilde, &Btilde[ j*k ], &gamma( 0,j ), ldC );
  }
}

void LoopOne( int m, int n, int k, double *Atilde, double *MicroPanelB, double *C, int ldC )
{
  for ( int i=0; i<m; i+=MR ) {
    int ib = min( MR, m-i );
    Gemm_6x8Kernel_Packed( k, &Atilde[ i*k ], MicroPanelB, &gamma( i,0 ), ldC );
  }
}

void PackMicroPanelA_MRxKC( int m, int k, double *A, int ldA, double *Atilde )
/* Pack a micro-panel of A into buffer pointed to by Atilde. 
   This is an unoptimized implementation for general MR and KC. */
{
  /* March through A in column-major order, packing into Atilde as we go. */

  if ( m == MR )   /* Full row size micro-panel.*/
    for ( int p=0; p<k; p++ ) 
      for ( int i=0; i<MR; i++ )
	*Atilde++ = alpha( i, p );
  else /* Not a full row size micro-panel.  We pad with zeroes. */
    for ( int p=0; p<k; p++ ) {
      for ( int i=0; i<m; i++ )
	*Atilde++ = alpha( i, p );
      for ( int i=m; i<MR; i++ )
	*Atilde++ = 0.0;
    }
}

void PackBlockA_MCxKC( int m, int k, double *A, int ldA, double *Atilde )
/* Pack a MC x KC block of A.  MC is assumed to be a multiple of MR.  The block is 
   packed into Atilde a micro-panel at a time. If necessary, the last micro-panel 
   is padded with rows of zeroes. */
{
  for ( int i=0; i<m; i+= MR ){
    int ib = min( MR, m-i );
    PackMicroPanelA_MRxKC( ib, k, &alpha( i, 0 ), ldA, Atilde );
    Atilde += ib * k;
  }
}

void PackMicroPanelB_KCxNR( int k, int n, double *B, int ldB, double *Btilde )
/* Pack a micro-panel of B into buffer pointed to by Btilde. 
   This is an unoptimized implementation for general KC and NR. */
{
  /* March through B in row-major order, packing into Btilde as we go. */

  if ( n == NR ) /* Full column width micro-panel.*/
    for ( int p=0; p<k; p++ ) 
      for ( int j=0; j<NR; j++ )
	*Btilde++ = beta( p, j );
  else /* Not a full row size micro-panel.  We pad with zeroes. */
    for ( int p=0; p<k; p++ ) {
      for ( int j=0; j<n; j++ )
	*Btilde++ = beta( p, j );
      for ( int j=n; j<NR; j++ )
	*Btilde++ = 0.0;
    }
}

void PackPanelB_KCxNC( int k, int n, double *B, int ldB, double *Btilde )
/* Pack a KC x NC panel of B.  NC is assumed to be a multiple of NR.  The block is 
   packed into Btilde a micro-panel at a time. If necessary, the last micro-panel 
   is padded with columns of zeroes. */
{
  for ( int j=0; j<n; j+= NR ){
    int jb = min( NR, n-j );
    
    PackMicroPanelB_KCxNR( k, jb, &beta( 0, j ), ldB, Btilde );
    Btilde += k * jb;
  }
}


#include <immintrin.h>

void Gemm_6x8Kernel_Packed( int k,
		        double *BlockA, double *PanelB, double *C, int ldC )
{
  __m256d gamma_0_0123 = _mm256_loadu_pd( &gamma( 0,0 ) );
  __m256d gamma_0_4567 = _mm256_loadu_pd( &gamma( 0,4 ) );
  __m256d gamma_1_0123 = _mm256_loadu_pd( &gamma( 1,0 ) );
  __m256d gamma_1_4567 = _mm256_loadu_pd( &gamma( 1,4 ) );
  __m256d gamma_2_0123 = _mm256_loadu_pd( &gamma( 2,0 ) );
  __m256d gamma_2_4567 = _mm256_loadu_pd( &gamma( 2,4 ) );
  __m256d gamma_3_0123 = _mm256_loadu_pd( &gamma( 3,0 ) );
  __m256d gamma_3_4567 = _mm256_loadu_pd( &gamma( 3,4 ) );
  __m256d gamma_4_0123 = _mm256_loadu_pd( &gamma( 4,0 ) );
  __m256d gamma_4_4567 = _mm256_loadu_pd( &gamma( 4,4 ) );
  __m256d gamma_5_0123 = _mm256_loadu_pd( &gamma( 5,0 ) );
  __m256d gamma_5_4567 = _mm256_loadu_pd( &gamma( 5,4 ) );

  __m256d alpha_i0_p;
  __m256d alpha_i1_p;
   	
  for ( int p=0; p<k; p++ ){
    /* load beta( p, 0:7 ) */
    __m256d beta_p_0123 = _mm256_loadu_pd( PanelB );
    __m256d beta_p_4567 = _mm256_loadu_pd( PanelB+4 );

    /* load alpha( 0, p ); update gamma( 0, 0:7 ) */
    alpha_i0_p = _mm256_broadcast_sd( BlockA );
    alpha_i1_p = _mm256_broadcast_sd( BlockA+1 );
    gamma_0_0123 = _mm256_fmadd_pd( alpha_i0_p, beta_p_0123, gamma_0_0123 );
    gamma_0_4567 = _mm256_fmadd_pd( alpha_i0_p, beta_p_4567, gamma_0_4567 );

    /* load alpha( 1, p ); update gamma( 1, 0:7 ) */
    gamma_1_0123 = _mm256_fmadd_pd( alpha_i1_p, beta_p_0123, gamma_1_0123 );
    gamma_1_4567 = _mm256_fmadd_pd( alpha_i1_p, beta_p_4567, gamma_1_4567 );

    /* load alpha( 2, p ); update gamma( 2, 0:7 ) */
    alpha_i0_p = _mm256_broadcast_sd( BlockA+2 );
    alpha_i1_p = _mm256_broadcast_sd( BlockA+3 );
    gamma_2_0123 = _mm256_fmadd_pd( alpha_i0_p, beta_p_0123, gamma_2_0123 );
    gamma_2_4567 = _mm256_fmadd_pd( alpha_i0_p, beta_p_4567, gamma_2_4567 );

    /* load alpha( 3, p ); update gamma( 3, 0:7 ) */
    gamma_3_0123 = _mm256_fmadd_pd( alpha_i1_p, beta_p_0123, gamma_3_0123 );
    gamma_3_4567 = _mm256_fmadd_pd( alpha_i1_p, beta_p_4567, gamma_3_4567 );

    /* load alpha( 4, p ); update gamma( 4, 0:7 ) */
    alpha_i0_p = _mm256_broadcast_sd( BlockA+4 );
    alpha_i1_p = _mm256_broadcast_sd( BlockA+5 );
    gamma_4_0123 = _mm256_fmadd_pd( alpha_i0_p, beta_p_0123, gamma_4_0123 );
    gamma_4_4567 = _mm256_fmadd_pd( alpha_i0_p, beta_p_4567, gamma_4_4567 );

    /* load alpha( 5, p ); update gamma( 5, 0:7 ) */
    gamma_5_0123 = _mm256_fmadd_pd( alpha_i1_p, beta_p_0123, gamma_5_0123 );
    gamma_5_4567 = _mm256_fmadd_pd( alpha_i1_p, beta_p_4567, gamma_5_4567 );
    
    BlockA += MR;
    PanelB += NR;
  }

  /* Store the updated results.  This should be done more carefully since
     there may be an incomplete micro-tile. */
  _mm256_storeu_pd( &gamma(0,0), gamma_0_0123 );
  _mm256_storeu_pd( &gamma(0,4), gamma_0_4567 );
  _mm256_storeu_pd( &gamma(1,0), gamma_1_0123 );
  _mm256_storeu_pd( &gamma(1,4), gamma_1_4567 );
  _mm256_storeu_pd( &gamma(2,0), gamma_2_0123 );
  _mm256_storeu_pd( &gamma(2,4), gamma_2_4567 );
  _mm256_storeu_pd( &gamma(3,0), gamma_3_0123 );
  _mm256_storeu_pd( &gamma(3,4), gamma_3_4567 );
  _mm256_storeu_pd( &gamma(4,0), gamma_4_0123 );
  _mm256_storeu_pd( &gamma(4,4), gamma_4_4567 );
  _mm256_storeu_pd( &gamma(5,0), gamma_5_0123 );
  _mm256_storeu_pd( &gamma(5,4), gamma_5_4567 );
}
