/*imports*/
#include "stdlib.h"
#include "string.h"
#include <math.h>
#include "conf.h"

/*definitions*/
#define PI 3.14159265358979323846264338327950288

/*
* Implementation of Hanning Window function in C
* HANNING(N) returns the N-point symmetric Hanning window in a column vector.
* 
* HANNING(N,'symmetric') returns the same result in HANNING(N)
* 
* HANNING(N,'periodic') returns the N-point periodic Hanning window, and includes the first 
* zero-weighted window sample
*
* itype = 1 --> periodic
* itype = 0 --> symmetric
*
*/

float *hanning(Int32 N, Int16 itype)
{
    Int32 half, i, idx, n;
    float *w;

    w = (float *)calloc(N, sizeof(float));
    memset(w, 0, N * sizeof(float));

    if (itype == 1) //periodic function
        n = N - 1;
    else
        n = N;

    if (n % 2 == 0)
    {
        half = n / 2;
        for (i = 0; i < half; i++) //CALC_HANNING   Calculates Hanning window samples.
            w[i] = 0.5 * (1 - cos(2 * PI * (i + 1) / (n + 1)));

        idx = half - 1;
        for (i = half; i < n; i++)
        {
            w[i] = w[idx];
            idx--;
        }
    }
    else
    {
        half = (n + 1) / 2;
        for (i = 0; i < half; i++) //CALC_HANNING   Calculates Hanning window samples.
            w[i] = 0.5 * (1 - cos(2 * PI * (i + 1) / (n + 1)));

        idx = half - 2;
        for (i = half; i < n; i++)
        {
            w[i] = w[idx];
            idx--;
        }
    }

    if (itype == 1)
    { //periodic function
        for (i = N - 1; i >= 1; i--)
            w[i] = w[i - 1];
        w[0] = 0.0;
    }
    return (w);
}

/*
* Implementation of Zeros function in C
* 
* X = zeros(4) --> create 4x4 matrix of zeros
* X = zeros(3,2) --> create 3x2 matrix of zeros
*
* So for general purposes,
* zeros(i,j) --> i=row size, j=column size
*/

float **zeros(Int32 rowSize, Int32 columnSize)
{
    Int32 i, j;
    float **arr = (float **)malloc(rowSize * sizeof(float *));
    for (i = 0; i < rowSize; i++)
    {
        arr[i] = (float *)malloc(columnSize * sizeof(float));
    }
    for (i = 0; i < rowSize; i++)
    {
        for (j = 0; j < columnSize; j++)
        {
            arr[i][j] = 0;
        }
    }
    return arr;
}

//one dimensional ones array for simplicity
float *zeros1d(Int32 length){
    Int32 i;
    float *arr;
    arr = (float*)calloc(length*sizeof(float));
    memset(arr, 0, length * sizeof(float));
    for(i=0;i<length;i++){
        arr[i]=0;
    }
    return arr;
}

/*
* Implementation of Ones function in C
* 
* X = ones(4) --> create 4x4 matrix of ones
* X = ones(3,2) --> create 3x2 matrix of ones
*
* So for general purposes,
* ones(i,j) --> i=row size, j=column size
*/

float **ones(Int32 rowSize, Int32 columnSize)
{
    Int32 i, j;
    float **arr = (float **)malloc(rowSize * sizeof(float *));
    for (i = 0; i < rowSize; i++)
    {
        arr[i] = (float *)malloc(columnSize * sizeof(float));
    }
    for (i = 0; i < rowSize; i++)
    {
        for (j = 0; j < columnSize; j++)
        {
            arr[i][j] = 1;
        }
    }
    return arr;
}

//one dimensional ones array for simplicity
float *ones1d(Int32 length){
    Int32 i;
    float *arr;
    arr = (float*)calloc(length*sizeof(float));
    memset(arr, 0, length * sizeof(float));
    for(i=0;i<length;i++){
        arr[i]=1;
    }
    return arr;
}

/*
* Implementation of nextpow2 function in MATLAB in C
* eg: parameter-->3 returns-->2 because 3<=4 && 4 = 2^2
* eg: parameter-->-3 returns-->2 because abs(-3)<=4 && 4 = 2^2
*/

Int32 nextpow2(Int32 n){
    return ceil(log2(abs(n)));
}

/*
* Implementation of min function
* minimum value between two numbers
*/

float min(float n1, float n2)
{
    return n1 <= n2 ? n1 : n2;
}

/*
* Implementation of max function
* maximum value between two numbers
*/

float max(float n1, float n2)
{
    return n1 >= n2 ? n1 : n2;
}

/*
* Implementation of flipud function of MATLAB in C
* flipud function flips rows in the up-down direction
* 
* Here, we are passing 2D array therefore we need to 
* specify the row size and column size
*
* We can find the row size and column size by following formulae
* int columnSize = sizeof((*arr)[0]);
* int rowSize = sizeof(*arr)/sizeof((*arr)[0]);
*/
float **flipud(int m, int n, int arr[][n]){
    Int32 i,j;
    Int32 columnSize = n;
    Int32 rowSize = m;

    float **flippedarr = (float **)malloc(rowSize * sizeof(float *));
    for (i = 0; i < rowSize; i++)
    {
        flippedarr[i] = (float *)malloc(columnSize * sizeof(float));
    }
    for (i = 0; i < rowSize; i++)
    {
        for (j = 0; j < columnSize; j++)
        {
            flippedarr[i][j] = arr[rowSize-i-1][j];
        }
    }
    return flippedarr;
}

/*
* Implementation of FFT
* FFT - Fast Fourier Transform
*/
void fft(complex *X, int M)
         /* X is an array of N = 2**M complex points. */
	   /* M = log2(N), N is the number of points */
  {
  complex temp1;          /* temporary storage complex variable */
  complex W;              /* e**(-j 2 PI/ N)                    */
  complex U;              /* Twiddle factor W**k                */
  int i,j,k;              /* loop indexes                       */
  int id;                 /* Index of lower point in butterfly  */
  int N = 1 << M;         /* Number of points for FFT           */
  int N2 = N/2;
  int L;                  /* FFT stage                          */
  int LE;           /* Number of points in sub DFT at stage L,  */
                    /* and offset to next DFT in stage          */
  int LE1;                /* Number of butterflys in one DFT at */
                          /*  stage L. Also is offset to lower  */
                          /*  point in butterfly at stage L     */

/*  Rearrange input array in bit-reversed order                 */
/*                                                              */
/*    The index j is the bit reversed value of i.  Since 0 -> 0 */
/*  and N-1 -> N-1 under bit-reversal, these two reversals are  */
/*  skipped.                                                    */

  j = 0;
  for(i=1; i<(N-1); i++)
    {
//  Increment bit-reversed counter for j by adding 1 to msb and */
//   propagating carries from left to right.                    */
                                            
    k = N2;               // k is 1 in msb, 0 elsewhere         

//  Propagate carry from left to right                          

    while(k<=j)           // Propagate carry if bit is 1        
      {
      j = j - k;          // Bit tested is 1, so clear it.      
      k = k/2;            // Set up 1 for next bit to right.    
      }
    j = j+k;              // Change 1st 0 from left to 1        
//  Swap samples at locations i and j if not previously swapped.

    if(i<j)               // Test if previously swapped.
      {
      temp1.real = (X[j]).real;
      temp1.imag = (X[j]).imag;
      (X[j]).real = (X[i]).real;
      (X[j]).imag = (X[i]).imag;
      (X[i]).real = temp1.real;
      (X[i]).imag = temp1.imag;
      }
    }

/*==============================================================*/
//  Do M stages of butterflys

  for(L=1; L<= M; L++)
    {
    LE = 1 << L;            //  LE = 2**L = points is sub DFT 
    LE1 = LE/2;             // Number of butterflys in sub-DFT
    U.real = 1.0;
    U.imag = 0.0;           // U = 1 + j 0                     
    W.real = cos(PI/LE1);
    W.imag = - sin(PI/LE1); // W = e**(-j 2 PI/LE)
//  Do butterflys for L-th stage

    for(j=0; j<LE1; j++)    //Do the LE1 butterflys per sub DFT
      {
//Compute butterflys that use same W**k
      for(i=j; i<N; i += LE)
        {
        id = i + LE1;      // Index of lower point in butterfly
        temp1.real = (X[id]).real*U.real - (X[id]).imag*U.imag;
        temp1.imag = (X[id]).imag*U.real + (X[id]).real*U.imag;

        (X[id]).real = (X[i]).real - temp1.real;
        (X[id]).imag = (X[i]).imag - temp1.imag;

        (X[i]).real  = (X[i]).real + temp1.real;
        (X[i]).imag  = (X[i]).imag + temp1.imag;
        }
    //Recursively compute W**k as W*W**(k-1) = W*U
      temp1.real = U.real*W.real - U.imag*W.imag;
      U.imag = U.real*W.imag + U.imag*W.real;
      U.real = temp1.real;
      }
    }
  return;
  }
