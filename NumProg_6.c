#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/times.h>
#include <limits.h>

/*1910049 Ishizeki Kota*/

#ifndef SIZE
#define SIZE	(1024)
#endif
#ifndef LIMIT
#define LIMIT	(200)
#endif


double	matA[SIZE][SIZE], vecB[SIZE], vecX[SIZE], vecT[SIZE];
int band[SIZE][2];  /*[0]には、左からいくつ0が続き、[1]には、右からいくつ0が続くかを記録する*/
int	n, verbose;

void Simal_Read ( double a[SIZE][SIZE], double b[SIZE] )	/* 連立方程式の読み込み */
{
  int	i, j;

  for ( i=0; i<n; i++ ) {
    for ( j=0; j<n; j++ ) {
      scanf( "%lf", &(a[i][j]) );
    }
    scanf( "%lf", &(b[i]) );
  }
}

void Simal_Write ( double a[SIZE][SIZE], double b[SIZE] )	/* 連立方程式の書き出し */
{
  int	i, j;

  for ( i=0; i<n; i++ ) {
    for ( j=0; j<n; j++ ) {
      printf( "%11.3le ", a[i][j] );
    }
    printf( "  %11.3le\n", b[i] );
  }
}

void V_Write ( double a[SIZE] )		       	/* ベクトルの書き出し */
{
  int i;
  for (i=0; i<n; i++) {
    printf( " %15.8le", a[i] );
  }
  printf( "\n" );
}

int main( int argc, char **argv ){

  int		i, j, k, cont, cc, dig, map;
  double 	C, Sum, prec, p;
  struct tms	tfrom, tto;
  char		ch, *ptr;
  int     cal_times = 0;
   
  printf("1910049 Ishizeki Kota\n");

  n = 0;
  verbose = 0;
  map = 0;
  C = 0;
  prec = 1e-5;

  while (( ch = getopt(argc,argv,"vms:i:p:")) != -1 ) {
    if ( ch == 'v' ) { verbose = 1; }		/* おしゃべり */
    if ( ch == 'm' ) { map = 1; }		/* マップ表示選択 */
    if ( ch == 's' ) { n = atoi( optarg ); }	/* 連立元数 Size */
    if ( ch == 'i' ) { C = atof( optarg ); }	/* 初期値 */
    if ( ch == 'p' ) {				/* 要求精度 (桁数) */
      j = atoi( optarg );
      prec = 1.0;
      for (i=0; i<j; i++) {
	prec *= 10.0;
      }
      prec = 1.0 / prec;
    }
  }

  if ( n == 0 ) {
    printf( "Size: " );
    scanf( "%d", &n );
    printf( "\n" );
  }
  if ( n > SIZE ) {			/* エラー処理 */
    printf( "error: Size should be less than or equal to %d.\n", SIZE );
    return 1;
  }

  for (i=0; i<n; i++) {
    vecX[i] = C;			/* 初期値代入 */
  }

  Simal_Read( matA, vecB );		/* 問題の読み込み */
  if (verbose) {
    Simal_Write( matA, vecB );		/* 清書 */
  }

  /*matAについて非ゼロ要素の範囲を見つける*/
  /*j行目についてband帯の大きさを求める*/
  for (j = 0; j < n; j++)
  {
    //左から調べる
    for(k = 0; k < n; k++){   //非ゼロ要素が左から詰まっているのなら0が入る。0がある数だけbandに入っていく
      if (matA[j][k] != 0){  //左から非ゼロ要素を探していき、ゼロになったところでループを抜ける
        break;
      }
      band[j][0]++;
    }

    for(k = n-1; k > 0; k--){   //非ゼロ要素が左から詰まっているのなら0が入る。0がある数だけbandに入っていく
      if (matA[j][k] != 0){
        break;
      }
      band[j][1]++;
    }

  }
  
  /*非ゼロ要素について範囲を絞れているか確認する*/
  /*
  printf("about band\n");
  for(j = 0; j < n; j++){
    printf("band[%d][0]:%d\n", j, band[j][0]);
  }
  for(j = 0; j < n; j++){
    printf("band[%d][1]:%d\n", j, band[j][1]);
  }
  */


  for (j=0; j<n; j++) {			/* 係数マトリクスの下準備 */
    C = matA[j][j];
    for (k=0; k<n; k++) {
      if (j != k) {
	      matA[j][k] = -matA[j][k]/C;
      }
    }
    vecB[j] = vecB[j]/C;
  }

  times( &tfrom );			/* 計算時間計測開始 */

  cc = 0;
  cont = 1;
  while ( cont && (cc<LIMIT) ) {
    if (!map) {
      printf( "%03d: ", cc );
      V_Write( vecX );			/* 暫定解の書き出し */
    }


    /*ここの計算に要するj,kについて限定してやる。あらかじめ0のある部分はループから除く*/
    for (j=0; j<n; j++) {
      Sum = vecB[j];

      for (k = band[j][0]; k < (n - band[j][1]); k++) {   //ループを回す範囲を限定する
	      if (j != k) {
	        Sum += matA[j][k] * vecX[k];	/* 積和計算 */
          cal_times++;
	      }
      }
      vecT[j] = Sum;
    }
    cc++; //修正計算を行うたびに1増える、解の有効桁数

    if (map) {
      printf( "%03d: ", cc );
    }
    cont = 0;
    for (i=0; i<n; i++) {		/* 収束状況の確認 */
      p = 999;
      if ( vecT[i]!=0 && (p=fabs((vecT[i]-vecX[i])/vecT[i])) > prec ) {
	      cont = 1; 
      }
      if (map) {			/* 有効桁数マップの作成 */
	      dig = (long)( -log(p)/log(10.0)+0.5 );
	      if (dig > 35) { dig = 35; }
	      if (dig < 0) { dig = 0; }
	      printf( "%c", (dig<=0) ? '.' : ((dig<10)?'0'+dig:'A'+dig-10) );
	      if (i%100 == 99) {
	        printf("\n     ");
	      }
      }
      vecX[i] = vecT[i];
    }
    if (map) {
      printf("\n");
    }
  }

  times( &tto );			/* 計算時間計測終了 */
  printf( "%03d: ", cc );
  V_Write( vecX );			/* 最終解の書き出し */
  printf("cal_times:%d\n", cal_times);
  printf( "    user cpu time: %.2lf [sec]\n\n",
	  (double)(tto.tms_utime-tfrom.tms_utime)/sysconf(_SC_CLK_TCK) );

  return 0;
}