#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

/* 1910049 Ishizeki Kota */

int verbose=0, query=1;

double func_R( double t){
	//return 10*sin(sqrt(0.1)*t + (M_PI/2) );
	return sin(2*t + (M_PI / 2));
}

double func( double t, double x, double v )
{
  //return -t*v - 20*x;
  return - 4.0 * x;
}

double vtox( double t, double x, double v )
{
  return v;
}


int main (int argc, char **argv) {
  double  t, tfrom, tto, h;
  double  x, a, b, c, d, delta;
  double  v, va, vb, vc, vd, vdelta;
  int	  n, interval=1;
  char    ch;

  printf("1910049 Ishizeki Kota\n\n");

  while (( ch = getopt(argc,argv,"vp:i:")) != -1 ) {
    if ( ch == 'v' ) {  verbose = 1;  }
    if ( ch == 'p' ) {  query = 0;
      sscanf( optarg, "%lf,%lf,%lf,%lf,%lf", &tfrom, &x, &v, &tto, &h );  }
    if ( ch == 'i' ) {  interval = atol( optarg );  }
  }

  if ( query ) {
    printf("t0 x0 v0 tto tinterval> ");
    scanf( "%lf %lf %lf %lf %lf", &tfrom, &x, &v, &tto, &h );
  }

  n = 0;
  t = tfrom;
  while ( t <= tto ) {
    a  = h * vtox( t,         x,         v          );
    va = h * func( t,         x,         v          );
    b  = h * vtox( t + h/2.0, x + a/2.0, v + va/2.0 );
    vb = h * func( t + h/2.0, x + a/2.0, v + va/2.0 );
    c  = h * vtox( t + h/2.0, x + b/2.0, v + vb/2.0 );
    vc = h * func( t + h/2.0, x + b/2.0, v + vb/2.0 );
    d  = h * vtox( t + h,     x + c,     v + vc     );
    vd = h * func( t + h,     x + c,     v + vc     );
    delta  = ( a  + 2.0*b  + 2.0*c  + d  ) / 6.0;
    vdelta = ( va + 2.0*vb + 2.0*vc + vd ) / 6.0;

    if ( n % interval == 0 ) {
      //printf( "%20.14lf %20.14lf %20.14lf %20.14lf %20.14lf\n", t, x, v, delta, vdelta );
      printf( "%20.14lf %20.14lf %20.14lf %20.14lf %20.14lf %20.14lf\n", t, x, func_R(t), v, delta, vdelta );
    }
    n++;
    t = tfrom + h * n;
    x += delta;
    v += vdelta;
  }

  return 0;
}