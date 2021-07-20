#include <math.h>
#include "stdlib.h"
#include "ezdsp5535_aic3204_dma.h"
#include "functions.h"

float *frequency_shaper(float *x, Int16 fs, float h_th, float p_th, float Lower_limit, float Upper_limit)
{
    /*
    * x - an input signal
    * fs - the sampling frequency of the input signal
    * h_th - threshold of hearing
    * p_th - threshold of pain
    * Lower_limit - Lower limit of the hearing loss range
    * Upper_limit - Upper limit of the hearing loss range
    */

    //iterator
    Int16 i = 0;

    float first = Lower_limit / 2;
    float second = Lower_limit;
    float third = Upper_limit;
    float fourth = Upper_limit + (fs / 2 - Upper_limit) / 2;

    Int16 x_length = sizeof(x) / sizeof(x[0]);
    Int16 n = nextpow2(x_length);
    Int32 N = (Int32)(pow(2, n));
    float T = 1 / (float)fs;

    complex *X = (complex *)malloc(N * sizeof(complex));
    complex *Y_shaped = (complex *)malloc(N * sizeof(complex));
    float *gain = zeros1d(N);

    //float firstC = (.3 * (h_th - 1)) / first;
    Int16 k = 0;

    float secondC = db2mag(1) + exp(5 / 2);
    float secondC2 = (float)(second - first) / 5;

    //float thirdC2 = (third - second) / 5;

    float fourthC = h_th;
    float fourthC2 = (float)(fourth - third) / 5;

    float fifthC = (float)(3 * h_th) / 4 + 1;
    float fifthC2 = (float)((float)fs / 2 - fourth) / 5;

    /*temporary variable*/
    complex *tmp;
    complex* Y = (complex*)malloc(N*sizeof(complex));
    complex* y_shaped = (complex*)malloc(N*sizeof(complex));
    float* y_shaped_clipped = (float*)malloc(x_length*sizeof(float));
    float* y_clipped = (float*)malloc(x_length*sizeof(float));
    
    fft(X, N);
    
    for (i = 0; i < N; i++)
    {
        Y_shaped[i].real = 0;
        Y_shaped[i].imag = 0;
    }

    //sets gain for the first stage of frequencies
    while ((float)k / N <= (float)first / fs)
    {
        gain[k] = db2mag(1) + exp(((float)(k - 1) / (N * T)) - (float)first / 2) / (float)((float)first / 5);
        gain[N - k - 1] = gain[k];
        Y_shaped[k].real = X[k].real * gain[k];
        Y_shaped[k].imag = X[k].imag * gain[k];
        Y_shaped[N - k - 1].real = X[N - k - 1].real * gain[N - k - 1];
        Y_shaped[N - k - 1].imag = X[N - k - 1].imag * gain[N - k - 1];
        if (sqrt(pow(Y_shaped[k].real, 2) + pow(Y_shaped[k].imag, 2)) > p_th)
        {
            Y_shaped[k].real = p_th;
            Y_shaped[k].imag = 0;
            Y_shaped[N - k - 1].real = p_th;
            Y_shaped[N - k - 1].imag = 0;
        }
        k++;
    }

    //sets the gain for the second stage of frequencies
    while ((float)k / N <= (float)second / fs)
    {
        gain[k] = h_th + (secondC - h_th) * exp(-(float)((float)(k - 1) / (N * T) - first) / secondC2);
        gain[N - k - 1] = gain[k];
        Y_shaped[k].real = X[k].real * gain[k];
        Y_shaped[k].imag = X[k].imag * gain[k];
        if (sqrt(pow(Y_shaped[k].real, 2) + pow(Y_shaped[i].imag, 2)) > p_th)
        {
            Y_shaped[k + 1].real = p_th;
            Y_shaped[k + 1].imag = 0;
            Y_shaped[N - k - 1].real = p_th;
            Y_shaped[N - k - 1].imag = 0;
        }
        k++;
    }

    //sets the gain for the third stage of frequencies
    while ((float)k / N <= (float)third / fs)
    {
        gain[k] = h_th;
        gain[N - k - 1] = gain[k];
        Y_shaped[k].real = X[k].real * gain[k];
        Y_shaped[k].imag = X[k].imag * gain[k];
        if (sqrt(pow(Y_shaped[k].real, 2) + pow(Y_shaped[k].imag, 2)) > p_th)
        {
            Y_shaped[k].real = p_th;
            Y_shaped[k].imag = 0;
            Y_shaped[N - k - 1].real = p_th;
            Y_shaped[N - k - 1].imag = 0;
        }
        k++;
    }

    //sets the gain for the fourth stage of frequencies
    while ((float)k / N <= (float)fourth / fs)
    {
        gain[k] = h_th - ((float)fourthC / 4 - 1) * exp(((float)k / (N * T) - fourth) / fourthC2);
        gain[N - k - 1] = gain[k];
        Y_shaped[k].real = X[k].real * gain[k];
        Y_shaped[k].imag = X[k].imag * gain[k];
        Y_shaped[N - k - 1].real = X[N - k - 1].real * gain[N - k - 1];
        Y_shaped[N - k - 1].imag = X[N - k - 1].real * gain[N - k - 1];
        if (sqrt(pow(Y_shaped[k].real, 2) + pow(Y_shaped[k].imag, 2)) > p_th)
        {
            Y_shaped[k].real = p_th;
            Y_shaped[k].imag = 0;
            Y_shaped[N - k - 1].real = p_th;
            Y_shaped[N - k - 1].imag = 0;
        }
        k++;
    }

    //sets the gain for the fifth stage stage of frequencies
    while((float)k/N<=.5){
        gain[k] = db2mag(1) + (fifthC-1)*exp(-(float)((float)k/(N*T)-fourth)/fifthC2);
        gain[N-k-1] = gain[k];
        Y_shaped[k].real = X[k].real*gain[k];
        Y_shaped[k].imag = X[k].imag*gain[k];
        Y_shaped[N-k-1].real = X[N-k-1].real*gain[N-k-1];
        Y_shaped[N-k-1].imag = X[N-k-1].real*gain[N-k-1];
        if (sqrt(pow(Y_shaped[k].real, 2) + pow(Y_shaped[k].imag, 2)) > p_th)
        {
            Y_shaped[k].real = p_th;
            Y_shaped[k].imag = 0;
            Y_shaped[N - k - 1].real = p_th;
            Y_shaped[N - k - 1].imag = 0;
        }
        k++;
    }
    for(i=0;i<N;i++){
        Y[i].real = X[i].real * gain[i];
        Y[i].imag = X[i].imag * gain[i];
    }
    ifft(Y,N,tmp);
    for(i=0;i<N;i++){
        y_shaped[i].real = Y_shaped[i].real;
        y_shaped[i].imag = Y_shaped[i].imag;
    }
    ifft(y_shaped,N,tmp);
    for(i=0;i<x_length;i++){
        y_clipped[i] = Y[i].real;
        y_shaped_clipped[i] = y_shaped[i].real;
    }
    return y_shaped_clipped;
}