#include "stdlib.h"
#include "ezdsp5535_aic3204_dma.h"
#include <math.h>

void VAD(Int32 *input, Int16 fs, Int16 frameCounter, Int16 n_frame, Int16 *iVad, Int16 *fVad,
         Int16 *RMSEArray, Int16 *periodicityArray, Int16 *ratioArray, Int16 *thresholdArray,
         float energyMin, float energyMax, float deltaEmin,
         float deltaEmax, Int16 lambda, Int16 threshold, Int16 initialEnergyMin)
{

    Int16 initialNoiseDuration = 50 * fs / 1000;

    //iterator
    Int16 i = 0;

    //window
    float *window;
    complex *theFFT;
    float *fftSq;

    //variables
    Int16 periodicity;
    Int16 fftLength;

    //Initialize Energy Variables
    Int16 lowerEnergyVoiceBand = 2000;
    Int16 upperEnergyVoiceBand = 4000;

    //Initialize periodicity and energy ratio thresholds
    Int16 periodicityThreshold = 1;
    Int16 energyRatioThreshold = 10;

    //compute the voiced energy band ratio
    float energyBelow;
    float energyAbove;
    float energyRatio;

    //get the root mean square energy of the input frame
    Int16 size = sizeof(input) / sizeof(input[0]);
    theFFT = (complex *)malloc(size * sizeof(complex));
    fftSq = (float *)malloc(size * sizeof(float));

    Int16 total = 0;
    float RMSE;
    for (i = 0; i < size; i++)
    {
        total += input[i] * input[i];
    }
    RMSE = sqrt(total / size);
    RMSEArray[frameCounter] = RMSE;
    if ((frameCounter * n_frame) < initialNoiseDuration)
    {
        if (RMSE > energyMax)
        {
            //set max
            energyMax = 5 * RMSE;
        }
        if (energyMin == 0)
        {
            //initialize Emin
            energyMin = RMSE;
            initialEnergyMin = energyMin;
        }
        else
        {
            if (RMSE > energyMin)
            {
                energyMin = RMSE;
            }
        }
        iVad[frameCounter] = 0;
        fVad[frameCounter] = 0;
    }
    else
    {
        //commence actual VAD processing
        //update running energy threshold estimates(Emax,Emin)
        if (RMSE > energyMax)
        {
            //we use a small delta to gradually decrease Emax to compensate for anomalous spikes in energy
            energyMax = RMSE;
            deltaEmax = 1;
        }
        else
        {
            deltaEmax = 0.999;
        }

        if (RMSE < energyMin)
        {
            if (RMSE == 0)
            {
                energyMin = initialEnergyMin;
            }
            else
            {
                energyMin = RMSE;
            }
            deltaEmin = 1;
        }
        else
        {
            //we use a small delta scaling factor to prevent complications, arising from energy dips(anomalies).
            //This keeps Emin rising at a gradual, minimal rate
            deltaEmin = deltaEmin * 1.001;
        }

        //threshold computations. Lambda is the non-linear dynamic coefficient used to compute the threshold
        //in a way that makes it resistant and independent of variations in background noise.
        lambda = 1 - (energyMin / energyMax);
        threshold = ((1 - lambda) * energyMax) + (lambda * energyMin);

        //keep track of threshold values
        thresholdArray[frameCounter] = threshold;

        //get the periodicity of the frame
        periodicity = Periodicity(input, fs);

        //add to periodicity array
        periodicityArray[frameCounter] = periodicity;

        //get the ratio of the frequencies above of the frequencies above and below 2kHz in the voice
        //band (0-4kHz). We will then use the ratio of these energies to make decisions around the voicing of the frame
        window = hamming(size);

        //fft computation and normalization
        fftLength = 2 ^ nextpow2(size);
        for (i = 0; i < size; i++)
        {
            theFFT[i].real = window[i] * input[i];
            theFFT[i].imag = 0;
        }
        fft(theFFT, fftLength);
        for (i = 0; i < size; i++)
        {
            fftSq[i] = (theFFT[i].real * theFFT[i].real + theFFT[i].imag * theFFT[i].imag) / (float)fftLength;
        }
        
        //compute the voiced energy band ratio
    }
}