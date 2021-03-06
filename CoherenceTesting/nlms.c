#include <math.h>
#include "stdlib.h"
#include "ezdsp5535_aic3204_dma.h"
#include "functions.h"

float *NLMS(Int32 *noisy, float *vad, float *w1, float *xx, float M)
{
    Int32 i, j;
    float tmp;
    float mu;
    float c = 0.01;
    float alpha = 0.09;
    float *noise = awgn(vad, 20); //20 SNR
    //Int32 size_vad = sizeof(vad) / sizeof(vad[0]);
    Int32 size_noise = sizeof(noise) / sizeof(noise[0]);
    float *d = (float *)malloc((size_noise) * sizeof(float));
    Int32 Ns = sizeof(d) / sizeof(d[0]);
    float *xxx = (float *)malloc(M * sizeof(float));
    float *y = (float *)malloc(M * sizeof(float));
    float *e = (float *)malloc(M * sizeof(float));

    for (i = 0; i < (size_noise); i++)
    {
        d[i] = noise[i] + vad[i];
    }
    for (i = 0; i < Ns; i++)
    {
        for (j = 0; j < (Int32)M - 1; j++)
        {
            xxx[j] = xx[j + 1];
        }
        xxx[(Int32)M - 1] = noise[i];
        for (j = 0; j < M; j++)
        {
            y[i] += w1[j] * xxx[j];
        }
        e[i] = d[i] - y[i];
        for (j = 0; j < M; j++)
        {
            tmp += xxx[j] * xxx[j];
        }
        mu = alpha / (c + tmp);
        for (j = 0; j < M; j++)
        {
            w1[j] = w1[j] + (mu * e[i] * xxx[j]);
        }
    }
    return y;
}
