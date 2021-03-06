#include <math.h>
#include "ezdsp5535_aic3204_dma.h"
#include "stdlib.h"
#include "functions.h"
#include "coherence.h"

#define SAMPLE_RATE 8000L
#define n_frame_coherence 640
#define n_frame_adaptive 320 //20ms processing from size

#pragma DATA_ALIGN(sampleBufferL, 4)
Int16 sampleBufferL[AUDIO_IO_SIZE];
#pragma DATA_ALIGN(sampleBufferR, 4)
Int16 sampleBufferR[AUDIO_IO_SIZE];
short tempBuff[128];

void main_function()
{
    /*Initialize parameters for VAD*/
    //not sure
    Int16 numFrames = ceil(AUDIO_IO_SIZE / n_frame_adaptive);

    /*array for coherence*/
    float *coherence_output = (float *)malloc(AUDIO_IO_SIZE * sizeof(float));
    Int16 coherence_output_clipper = 5 * n_frame_adaptive / 4 - +(n_frame_adaptive / 4);
    float *coherence_output_new = (float *)malloc(coherence_output_clipper * sizeof(float));

    /*array for NLMS*/
    float *y1 = (float *)malloc(coherence_output_clipper * sizeof(float));

    /*frequency shaper*/
    float *freq_shaper;

    /*Initialize Energy Variables*/
    float energyMin = 0;
    float energyMax = 0;
    float deltaEmin = 1;
    float deltaEmax = 1;
    float lambda = 1;
    float threshold = 0;
    float initialEnergyMin = 0;

    /*Initialize parameters for NLMS adaptive filter*/
    float M = 70;
    float *xx = zeros1d(M);
    float *w1 = zeros1d(M);

    /*Initializing parameters for frequency shaper*/
    float H_th = db2mag(20);  //threshold of hearing
    float P_th = db2mag(70);  //threshold of pain
    Int16 Lower_limit = 2000; //lower limit of hearing loss range
    Int16 Upper_limit = 5000; //upper limit of hearing loss range

    /*Iterators*/
    Int16 i;

    EZDSP5535_init();

    EZDSP5535_SAR_init();

    printf("\n Getting audio signals \n");

    aic3204_hardware_init();

    aic3204_init();

    aic3204_dma_init();

    set_sampling_frequency_and_gain(SAMPLE_RATE, 0);

    while (1)
    {
        aic3204_read_block(sampleBufferL, sampleBufferR);

        coherence_output = COH_RealImag(sampleBufferL, sampleBufferR, SAMPLE_RATE); //not sure that the sample rate the one should come
        //we are going to have a continuous output, therefore
        for (i = 1 + (n_frame_adaptive / 4); i < 5 * n_frame_adaptive / 4; i++)
        {
            coherence_output_new[i - (1 + (n_frame_adaptive / 4))] = coherence_output[i];
        }
        //NLMS Adaptive Filter
        y1 = NLMS(coherence_output_new, w1, xx, M); //pass all as pointers
                                                    //Frequency Shaper
        //don't know the dimensions
        freq_shaper = freqencyshaper(y1, SAMPLE_RATE, H_th, P_th, Lower_limit, Upper_limit);
        aic3204_write_block(sampleBufferR, sampleBufferR);
    }
}