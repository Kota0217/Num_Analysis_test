#include <stdio.h>
#include <math.h>

/* 1910049 Ishizeki Kota */

#ifndef  SIZE
#define  SIZE	(16)
#endif


double matA[SIZE][SIZE], matE[SIZE][SIZE], vecX[SIZE];

void dump ( int n )
{
  int    j, k;
  for ( k=0; k<n; k++ ){
    printf( "\n" );
    for ( j=0; j<n; j++ ) printf( " %10.3le", matA[k][j] );
    printf(                    "    %10.3le   ", vecX[k] );
    for ( j=0; j<n; j++ ) printf( " %10.3le", matE[k][j] );
  }
  printf( "\n" );
}


int main( int argc, char **argv )
{
  double C;
  int	 n, i, j, k;
  int    ch, verbose;
  
  printf("1910049 Ishizeki Kota\n");

  verbose = 0;
  while (( ch = getopt(argc,argv,"vh")) != -1 ) {
    if ( ch == 'v' ) verbose = 1;
    if ( ch == 'h' ) {
      printf( "Options:\n" );
      printf( "\t-v: dumps intermediate results.\n" );
      printf( "\t-h: shows help message, as you see now.\n" );
      printf( "\n" );
      return 1;
    }
  }


/* データ入力 */
  scanf( "%d", &n );
  if ( n > SIZE ) {
    printf("N should be less than or equal to %d.\n", SIZE);
    return 1;
  }

  for ( i=0; i<n; i++ ) {
    for ( j=0; j<n; j++ ) {
      scanf( "%lf", &(matA[i][j]) ); }
    scanf( "%lf", &(vecX[i]) );
  }

  printf("Inputs:\n");
  
  for ( i=0; i<n; i++ ){
    for ( j=0; j<n; j++ ){
      printf( " %12.5le", matA[i][j] ); }
    printf( "  |  %12.5le\n", vecX[i] );
  }


/* 単位行列の設定 */
  for ( i=0; i<n; i++ ){
    for ( j=0; j<n; j++ ){ matE[i][j] = 0.0; }
    matE[i][i] = 1.0;
  }

/* Solve */
  if ( verbose ) dump( n );

  for ( i=0; i<n; i++ ){
    C = matA[i][i];

    //pivotが0、または0に近い行列であるときにそれを検出する
    if ( C == 0 ) {
      printf( "Ill conditioned matrix, A[%d][%d] = %le\n\n", i, i, C );
      return 1;
    }
    if (fabs(C) < 1e-4){
        printf("Ill conditioned matrix, A[%d][%d] = %le\n\n", i, i, C);
        return 1;
    }

    /*i 列目前半の処理*/
    for ( j=0; j<n; j++){
		matA[i][j] = matA[i][j] / C;
	}
	vecX[i] = vecX[i] / C;
		
	for (j=0; j<n; j++){
		matE[i][j] = matE[i][j] / C;
	}

    if ( verbose ) dump( n );

    for ( k=0; k<n; k++ ){
      if ( i != k ){
        C = matA[k][i];
	/********************************************/
	/**            i 列目後半の処理            **/
	/********************************************/
				for ( j=0; j<n; j++){
					matA[k][j] = matA[k][j] - C * matA[i][j];

				}
				vecX[k] = vecX[k] - C * vecX[i];
				for (j=0; j<n; j++){
					matE[k][j] = matE[k][j] - C * matE[i][j];
				}
      }
    }

    if ( verbose ) dump( n );
  }


/* 計算結果出力 */
  printf( "\nResults:\n" );
  for ( i=0; i<n; i++ ) printf( " %12.5le", vecX[i] );
  printf( "\nInverse:\n" );
  for ( i=0; i<n; i++ ){
    for ( j=0; j<n; j++ ) printf( " %12.5le", matE[i][j] );
    printf("\n");
  }
  printf("\n");

}