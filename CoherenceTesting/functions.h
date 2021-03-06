#ifndef FUNCTIONS_H_
#define FUNCTIONS_H_

#include "ezdsp5535_aic3204_dma.h"
#include "stdlib.h"
#include "string.h"
/****************************************************************
*				      Function Prototypes
****************************************************************/
float* hanning(Int32 N, Int16 itype);
float **zeros(Int32 rowSize, Int32 columnSize);
float *zeros1d(Int32 length);
float **ones(Int32 rowSize, Int32 columnSize);
float *ones1d(Int32 length);
Int32 nextpow2(Int32 n);
float min(float n1,float n2);
float max(float n1,float n2);
void fft(complex* X, int M);
float db2mag(float r);
void ifft( complex *v, int n, complex *tmp );
float *awgn(float *signal, float SNR);
#endif /*FUNCTIONS_H_*/
