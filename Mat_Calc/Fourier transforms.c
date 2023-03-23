#include<stdio.h> 
#include <complex.h> //complex.h を include する。 
#include <math.h>
	
#define N 101
#define K 60
#define M 1

/*プロトタイプ宣言*/
/*逆行列*/
void INV(int KK, double complex x[KK][KK]);
/*複素共役転置行列*/
void Conjtrans(int NN, int KK, double complex x[NN][KK], double complex y[NN][KK]);
/*行列の積*/
void PRD(int KK, int NN, int MM, double complex x[KK][NN], double complex y[NN][MM], double complex z[KK][MM]);

int main(){
	
	int i,j;
	double data[101][3];		//データファイル読み込み用
	char   p;					//csvファイルカンマ記号の一時格納用
	double re_A, im_A;			//実部，虚部一時格納用
	double complex Y[N][M];		//データ格納用
	double complex Y2[N][M];	//フィッティングデータ格納用
	double complex X1[N][K];	//行列X（基底行列）格納用
	double complex X2[K][N];	//途中計算結果格納用
	double complex X3[K][K];	//途中計算結果格納用
	double complex X4[K][N];	//途中計算結果格納用
	double complex a[K][M];		//係数ベクトル格納用
	FILE *fp; 	//読み込み用
	FILE *fp2;  //書き出し用

  FILE *gid = NULL;   //gnuplot用

	
/*外部データ読み込み部*/
	fp=fopen("data3.csv","r");

	printf("data\n");
	for(i=0;i<101;i++){
		for(j=0;j<3;j++){
			if(j!=2)fscanf(fp,"%lf,",&data[i][j]);
			if(j==2)fscanf(fp,"%lf",&data[i][j]);
			printf("%lf\t",data[i][j]);
		}
		fscanf(fp,"\n");
		printf("\n");
	}
	printf("\n\n");
	
/*ベクトルY作成*/
	printf("Y\n");
	for(i=0;i<N;i++){
		for(j=0;j<M;j++){
			Y[i][0]= data[i][1] + data[i][2]*I;
			printf("%f+%fi\n", creal(Y[i][j]), cimag(Y[i][j]));
		}
		printf("\n");
	}
	printf("\n\n");
	
/*行列X（基底行列）作成*/
	printf("X\n");
	for(i=0;i<N;i++){
		for(j=0;j<K;j++){
			re_A = cos(2*M_PI*1/N*(j-K/2)*i);
			im_A = sin(2*M_PI*1/N*(j-K/2)*i);
			X1[i][j]= re_A + im_A*I;
			printf("%lf+%lfi\t", creal(X1[i][j]), cimag(X1[i][j]));
		}
		printf("\n");
	}
	printf("\n\n");
	
/*一般化逆行列の計算*/
    //行列Xの転置
    Conjtrans(N,K,X1,X2);
    //行列X(転置)とXの積
    PRD(K,N,K,X2,X1,X3);
    //逆行列
  	INV(K,X3);
  	//逆行列と行列X(転置)の積
    PRD(K,K,N,X3,X2,X4);
    //ベクトルYとの積
    PRD(K,N,M,X4,Y,a);
    
/*フィッティングした波形の作成*/
	PRD(N,K,M,X1,a,Y2);
	
/*出力結果書き出し部*/
    fp2=fopen("data3_output.csv","w");
    for(i=0;i<N;i++){
    	fprintf(fp2,"%d,%lf,%lf\n",i,creal(Y2[i][0]),cimag(Y2[i][0]));
    }
    fclose(fp2);
/*以下，gnuplotによる描画部*/
  gid =popen("gnuplot","w");
  fprintf(gid, "set datafile separator ','\n");
  fprintf(gid," set xrange [0:100]\n");
  fprintf(gid," set yrange [-2:10]\n");
  fprintf(gid, "plot 'data3.csv' w p pt 6 ps 0.7 lt 1\n");
  fprintf(gid, "replot 'data3_output.csv' w p pt 2 ps 0.7 lt 1\n");
  fprintf(gid, "pause 10\n");
  fclose(gid);

}


/*行列と行列の積を求める関数*/
void PRD(int KK, int NN, int MM, double complex x[KK][NN], double complex y[NN][MM], double complex z[KK][MM]){
	int i,j,k;
	printf("PRD\n");
	for(i=0;i<KK;i++){
		for(j=0;j<MM;j++){
			z[i][j]=0;
			for (k=0;k<NN;k++) {
				z[i][j] += x[i][k] * y[k][j];
			}
		printf("%f+%fi\t", creal(z[i][j]), cimag(z[i][j]));
		}
		printf("\n");
	}
	printf("\n\n");
}

/*正方行列の逆行列を求める関数*/
void INV(int NN, double complex x[NN][NN]){
	int i,j,k;
	double complex c1[NN];
	double complex c;
	double complex II[NN][NN];
	
	printf("INV\n");
	
	for(i=0;i<NN;i++){
		for(j=0;j<NN;j++){
			II[i][j]=0 + 0*I;
			if(i==j) II[i][j]=1 + 0*I;
		}
	}
	for(i=0;i<NN;i++){
		for(j=0;j<NN;j++){
			printf("%lf+%lfi\t",creal(x[i][j]), cimag(x[i][j]));
		}
		printf("\n");
	}
	printf("\n\n");
	//（左下の成分を0にする）
	for(k=0;k<NN-1;k++){
		for(i=k+1;i<NN;i++){
			c=(x[i][k]/x[k][k]);
			for(j=k;j<NN;j++){
				x[i][j] = x[i][j]-x[k][j]*c;
			}
			for(j=0;j<NN;j++){
				II[i][j] = II[i][j]-II[k][j]*c;
			}
		}
	}

	//左下三角行列の確認
	for(i=0;i<NN;i++){
		for(j=0;j<NN;j++){
			printf("%lf+%lfi\t",creal(x[i][j]), cimag(x[i][j]));
		}
		printf("\n");
	}
	printf("\n\n");
	
	//（右上の成分を0にする）
	for(k=NN-1;k>0;k--){
		for(i=k-1;i>=0;i--){
			c=(x[i][k]/x[k][k]);
			for(j=k;j>=0;j--){
				x[i][j] = x[i][j]-x[k][j]*c;
			}
			for(j=0;j<NN;j++){
				II[i][j] = II[i][j]-II[k][j]*c;
			}
		}
	}
	//ここまで対角化完了
	
	//対角成分を1とする
	for(i=0;i<NN;i++){
				c1[i]=1/x[i][i];
				x[i][i]=c1[i]*x[i][i];
	}
	
	for(i=0;i<NN;i++){
		for(j=0;j<NN;j++){
			II[i][j]=c1[i]*II[i][j];
		}
	}

	//結果の確認（aが単位行列になっていることを確認する）
	for(i=0;i<NN;i++){
		for(j=0;j<NN;j++){
		printf("%lf+%lfi\t",creal(x[i][j]), cimag(x[i][j]));
		}
		printf("\n");
	}
	printf("\n\n");

	//Iが逆行列になっていることを確認する
	for(i=0;i<NN;i++){
		for(j=0;j<NN;j++){
			x[i][j]=II[i][j];
			printf("%lf+%lfi\t",creal(II[i][j]), cimag(II[i][j]));
		}
		printf("\n");
	}
	printf("\n\n");
}

//複素共役転置行列関数
void Conjtrans(int NN, int KK, double complex x[NN][KK], double complex y[KK][NN]){
		int i,j;
		for(i=0;i<NN;i++){
			for(j=0;j<KK;j++){
				y[j][i]=conj(x[i][j]);
			}
		}
		printf("TRP\n");
		for(i=0;i<KK;i++){
			for(j=0;j<NN;j++){
				printf("%lf+%lfi\t",creal(y[i][j]), cimag(y[i][j]));
			}
			printf("\n");
		}
		printf("\n\n");
}