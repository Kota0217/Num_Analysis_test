#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* 1910049 Ishizeki Kota */
/* シンプソン 3/8 公式 */

double fun( double x )
{
    //return x * exp(x);
    return tan(x)*tan(x);
}

double fun_R(double xfrom, double xto) /*解析的に解いた積分結果、xfromとxtoの代入値の差を返す*/
{
    /*
    xfrom = xfrom * exp(xfrom) - exp(xfrom);
    xto = xto * exp(xto) - exp(xto);
    */
    xfrom = tan(xfrom) - xfrom;
    xto = tan(xto) - xto;
    return xto - xfrom;
}


int main ( int argc, char **argv ){
    double  x, sum;
    double  xfrom, xto, hh;
    int	  xsteps, xidx;
    int i = 0;
    
    printf("1910049 Ishizeki Kota\n\n");

    printf("x: from to steps> ");
    if ( argc >= 4 ) {
        xfrom  = atof( argv[1] );
        xto    = atof( argv[2] );
        xsteps = atol( argv[3] );
        printf( "%lf %lf %d\n", xfrom, xto, xsteps );
    } else {
        scanf( "%lf %lf %d", &xfrom, &xto, &xsteps );
    }

    while (1)
    {
        if(i >= 1){
            printf("x: from to steps> ");
            scanf( "%lf %lf %d", &xfrom, &xto, &xsteps );
            printf( "%lf %lf %d\n", xfrom, xto, xsteps );
        }

        i++;


        if(xsteps <= 0){
            printf("xsteps <= 0\n");
            exit(1);
        }

        if( xsteps % 3 ) {
            while(xsteps % 3 != 0){
                printf( "Bad steps\n" );
                printf("ReEnter steps > ");
                scanf("%d", &xsteps);
                if (xsteps <= 0)
                {
                    printf("xsteps <= 0\n");
                    exit(1);
                }
            }
        }
  
        printf("Ans:\t%20.15f\n", fun_R(xfrom, xto));
        sum = 0;
        hh = ( xto - xfrom ) / xsteps;
        for (xidx=0; xidx<xsteps; xidx+=3) {
            x = (xfrom * (xsteps - xidx) + xto * xidx) / xsteps;	/* 内分点の計算 */
            sum += 2 * fun( x ) +  3 * ( fun( x + hh ) +  fun( x + 2 * hh ) );
        }
        sum += - fun( xfrom ) + fun( xto );
        sum *= 3.0 * (xto - xfrom) / ( 8.0 * xsteps );
        printf( "ans:\t%20.15lf\n", sum );
    }
}