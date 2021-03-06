#include "conf.h" /*Definition of complex variable structure*/
#include <math.h> /*Definitions of math library*/

//#define q 3 /*for 2^3 points*/
//#define N (1<<q) /*N-point FFT, iFFT*/
#define PI 3.14159265358979323846264338327950288

/* 
   fft(v,N):
   [0] If N==1 then return.
   [1] For k = 0 to N/2-1, let ve[k] = v[2*k]
   [2] Compute fft(ve, N/2);
   [3] For k = 0 to N/2-1, let vo[k] = v[2*k+1]
   [4] Compute fft(vo, N/2);
   [5] For m = 0 to N/2-1, do [6] through [9]
   [6]   Let w.re = cos(2*PI*m/N)
   [7]   Let w.im = -sin(2*PI*m/N)
   [8]   Let v[m] = ve[m] + w*vo[m]
   [9]   Let v[m+N/2] = ve[m] - w*vo[m]
 */
 
 void fft(complex *v, int n, complex *tmp){
 	if(n>1){ /*otherwise, do nothing and return*/
 		int k,m; 
 		complex z, w, *vo, *ve;
 		ve = tmp;
 		vo = tmp+n/2;
 		for(k=0; k<n/2;k++){
 			ve[k] = v[2*k];
 			vo[k] = v[2*k+1];
 		}
 		fft(ve, n/2, v); /*FFT on even-indexed elements of v[]*/
 		fft(vo, n/2, v); /*FFT on odd-indexed elements of v[]*/
 		for(m=0;m<n/2;m++){
 			w.real = cos(2*PI*m/(double)n); 
 			w.imag = -sin(2*PI*m/(double)n);
 			z.real = w.real*vo[m].real - w.imag*vo[m].imag; /*Re(w*vo[m])*/
 			z.imag = w.real*vo[m].imag + w.imag*vo[m].real; /*Im(w*vo[m])*/
 			v[  m  ].real = ve[m].real + z.real;
 			v[  m  ].imag = ve[m].imag + z.imag;
 			v[m+n/2].real = ve[m].real - z.real;
 			v[m+n/2].real = ve[m].imag - z.imag; 
 		}
 	}
 	return;
 }
 
 /* 
   ifft(v,N):
   [0] If N==1 then return.
   [1] For k = 0 to N/2-1, let ve[k] = v[2*k]
   [2] Compute ifft(ve, N/2);
   [3] For k = 0 to N/2-1, let vo[k] = v[2*k+1]
   [4] Compute ifft(vo, N/2);
   [5] For m = 0 to N/2-1, do [6] through [9]
   [6]   Let w.re = cos(2*PI*m/N)
   [7]   Let w.im = sin(2*PI*m/N)
   [8]   Let v[m] = ve[m] + w*vo[m]
   [9]   Let v[m+N/2] = ve[m] - w*vo[m]
 */
void
ifft( complex *v, int n, complex *tmp )
{
  if(n>1) {			/* otherwise, do nothing and return */
    int k,m;    
    complex z, w, *vo, *ve;
    ve = tmp; vo = tmp+n/2;
    for(k=0; k<n/2; k++) {
      ve[k] = v[2*k];
      vo[k] = v[2*k+1];
    }
    ifft( ve, n/2, v );		/* FFT on even-indexed elements of v[] */
    ifft( vo, n/2, v );		/* FFT on odd-indexed elements of v[] */
    for(m=0; m<n/2; m++) {
      w.real = cos(2*PI*m/(double)n);
      w.imag = sin(2*PI*m/(double)n);
      z.real = w.real*vo[m].real - w.imag*vo[m].imag;	/* Re(w*vo[m]) */
      z.imag = w.real*vo[m].imag + w.imag*vo[m].real;	/* Im(w*vo[m]) */
      v[  m  ].real = ve[m].real + z.real;
      v[  m  ].imag = ve[m].imag + z.imag;
      v[m+n/2].real = ve[m].real - z.real;
      v[m+n/2].imag = ve[m].imag - z.imag;
    }
  }
  return;
}
/**************************************************
 Print a vector of complexes as ordered pairs. 
**************************************************/
/*
static void
print_vector(
	     const char *title,
	     complex *x,
	     int n)
{
  int i;
  printf("%s (dim=%d):", title, n);
  for(i=0; i<n; i++ ) printf(" %5.2f,%5.2f ", x[i].Re,x[i].Im);
  putchar('\n');
  return;
}
*/

/*************************************************
 Implementation Example
**************************************************/
/*
int main(void)
{
  complex v[N], v1[N], scratch[N];
  int k;

  //Fill v[] with a function of known FFT:
  for(k=0; k<N; k++) {
    v[k].Re = 0.125*cos(2*PI*k/(double)N);
    v[k].Im = 0.125*sin(2*PI*k/(double)N);
    v1[k].Re =  0.3*cos(2*PI*k/(double)N);
    v1[k].Im = -0.3*sin(2*PI*k/(double)N);
  }
    
  // FFT, iFFT of v[]: 
  print_vector("Orig", v, N);
  fft( v, N, scratch );
  print_vector(" FFT", v, N);
  ifft( v, N, scratch );
  print_vector("iFFT", v, N);

  //FFT, iFFT of v1[]: 
  print_vector("Orig", v1, N);
  fft( v1, N, scratch );
  print_vector(" FFT", v1, N);
  ifft( v1, N, scratch );
  print_vector("iFFT", v1, N);

  exit(EXIT_SUCCESS);
}
*/