/*includes*/
#include "stdio.h"
#include "conf.h"
#include <math.h>

/*definitions*/
#define lambdaX 0.68 //forgetting factor for smoothing power spectrum
#define epsilon 10^(-12)
#define negImgConsFilter 0.05

/*
* Implementation of Coherence function
* x1,x2 --> Input Signals at two channels(vectors)
* fs --> sampling frequency
*/
void COH_RealImag(float *x1, float *x2, Int32 fs){
    //variable declaration
    Int32 i,j; //iteration values
    Int32 frameLength; //frame length
    Int32 frameShift; //frame shift
    Int32 FFT_LEN; //fft length
    Int32 lenS;//length of signal
    Int32 nFrame; //frame number --> initially 0
    Int32 iniFrameSample; //initial frame sample
    Int32 endFrameSample; //ending frame sample
    Int32 band1;
    Int32 rowSize, columnSize; //temporary variables to hold the array dimensions
    float *window; //hanning window --> In matlab we get a column vector here, we have an array which is a row vecotr
    float *P;
    float *limNeg;
    float *Frame1;
    float *Frame2;
    float *wFrame1;
    float *wFrame2;
    complex *X1;
    complex *X2;
    float *PX1X1;
    float *PX2X2;
    complex *PX1X2;
    complex *cohX;
    float *reCOH;
    float *imCOH;
    float *G1;
    float *G2;
    float *ind_neg;
    float *G;
    float *H;
    complex *HX1; //temporary variable
    float *enhSpeech_Frame_tmp;
    float *enhSpeech_Frame;
    float *enhanced_output;
    
    float* ptrFloat; //for freeing 2D arrays

    //find the frame length
    frameLength = floor(20*fs/1000);
    if(frameLength/(float)2!=floor(frameLength/(float)2)){
        frameLength = frameLength+1;
    }

    frameShift = floor(frameLength * 0.5);
    window = hanning(frameLength, 1); //we are going with symmetric
    FFT_LEN  = exp2(nextpow2(frameLength));

    //as x1,x2 are row vectors/arrays in here, we take the minimum value
    lenS = min((*(&x1+1)-x1),(*(&x2+1)-x2));
    nFrame = 0;
    iniFrameSample = 1;
    endFrameSample = iniFrameSample+frameLength-1;
    enhanced_output = (float*)malloc(lenS*sizeof(float));

    //Defining bands
    band1 = floor(1000*FFT_LEN/(float)fs); //1kHZ
    //Defining exponents of coherence function
    P = zeros1d(FFT_LEN/2);
    for(i=0;i<band1;i++){
        P[i] = 16; //since this is a column vector there is only one column
    }
    for(i=band1;i<FFT_LEN/2;i++){
        P[i] = 2;
    }

    //Defining thresholds for imaginary parts to consider negative noise
    limNeg = zeros1d(FFT_LEN/2);
    for(i=0;i<band1;i++){
        limNeg[i] = -0.1; //since this is a column vector there is only one column
    }
    for(i=band1;i<FFT_LEN/2;i++){
        limNeg[i] = -0.3;
    }

    //allocate memory for arrays
    Frame1 = (float*)malloc(frameLength, sizeof(float));
    Frame2 = (float*)malloc(frameLength, sizeof(float));

    wFrame1 = (float*)malloc(frameLength, sizeof(float));
    wFrame2 = (float*)malloc(frameLength, sizeof(float));

    X1 = (complex*)malloc(FFT_LEN, sizeof(complex));
    X2 = (complex*)malloc(FFT_LEN, sizeof(complex));

    PX1X1 = (float*)malloc(frameLength, sizeof(float));
    PX2X2 = (float*)malloc(frameLength, sizeof(float));
    PX1X2 = (complex*)malloc(frameLength, sizeof(complex));

    cohX = (complex*)malloc(frameLength, sizeof(complex));

    reCOH = (float*)malloc(FFT_LEN/2, sizeof(float));
    imCOH = (float*)malloc(FFT_LEN/2, sizeof(float));

    G1 = (float*)malloc(FFT_LEN/2,sizeof(float));
    G2 = ones1d(FFT_LEN/2);

    ind_neg = (float*)malloc(FFT_LEN/2,sizeof(float));

    G = (float*)malloc(FFT_LEN/2, sizeof(float));

    H = (float*)malloc(FFT_LEN * sizeof(float *));


    HX1 = (complex*)malloc(FFT_LEN*sizeof(complex));

    enhSpeech_Frame_tmp = (float*)malloc(FFT_LEN * sizeof(float));
    enhSpeech_Frame = (float*)malloc(frameLength*sizeof(float));

    while(endFrameSample<lenS){
        nFrame = nFrame + 1; //a new frame to process
        //get short-time magnitude and phase spectrum for each input channel
        for(i=iniFrameSample-1;i<endFrameSample;i++){
            Frame1[i] = x1[i];
            Frame2[i] = x2[i];
        }
        for(i=0;i<frameLength;i++){
            wFrame1[i] = Frame1[i] * window[i];
            wFrame2[i] = Frame2[i] * window[i];
        }
        for(i=0;i<frameLength;i++){
            X1[i].real = wFrame1[i];
            X1[i].imag = 0;
            X2[i].real = wFrame2[i];
            X2[i].imag = 0;
        }
        fft(X1, FFT_LEN);
        fft(X2, FFT_LEN);

        if(nFrame==1){
            for(i=0;i<frameLength;i++){
                PX1X1[i] = pow(X1[i].real,2)+pow(X1[i].imag,2);
                PX2X2[i] = pow(X2[i].real,2)+pow(X2[i].imag,2);
                PX1X2[i].real = X1[i].real*X2[i].real+X1[i].imag*X2[i].imag;
                PX1X2[i].imag = X1[i].real*X2[i].imag+X[i].imag*X2[i].real;
            }
        }else{
            for(i=0;i<frameLength;i++){
                PX1X1[i] = lambdaX*PX1X1[i]+(1-lambdaX)*pow(X1[i].real,2)+pow(X1[i].imag,2);
                PX2X2[i] = lambdaX*PX2X2[i]+(1-lambdaX)*pow(X2[i].real,2)+pow(X2[i].imag,2);
                PX1X2[i].real = lambdaX*PX1X2[i].real+(1-lambdaX)*X1[i].real*X2[i].real+X1[i].imag*X2[i].imag;
                PX1X2[i].imag = lambdaX*PX1X2[i].imag+(1-lambdaX)*X1[i].real*X2[i].imag+X[i].imag*X2[i].real;
            }
        }

        for(i=0;i<frameLength;i++){
            cohX[i].real = PX1X2[i].real/(float)(PX1X1[i]*PX2X2[i]);
            cohX[i].imag = PX1X2[i].imag/(float)(PX1X1[i]*PX2X2[i]);
        }

        for(i=1;i<FFT_LEN/2+1;i++){
            reCOH[i-1] = cohX[i].real;
            imCOH[i-1] = cohX[i].imag;
        }
        //for suppressing noise coming from angles about 90
        for(i=0;i<FFT_LEN/2;i++){
            G1[i] = 1 - abs(reCOH[i])*P[i];
        }
        //for suppressing noise coming from angles greater than 90 we have G2
        for(i=0;i<FFT_LEN/2;i++){
            ind_neg[i] = imCOH[i]<limNeg[i] ? 1:0;
        }
        for(i=0;i<FFT_LEN;i++){
            if(ind_neg[i]==1){
                G2[i]=negImgConsFilter;
            }
        }
        for(i=0;i<FFT_LEN/2;i++){
            G[i] = G1[i]*G2[i];
        }
        //fullband final filter
        for(i=0;i<FFT_LEN/2;i++){
            H[i] = G[i];
        }
        for(i=0;i<FFT_LEN/2;i++){
            H[i] = G[FFT_LEN/2-i-1];
        }
        //ifft and OLA
        for(i=0;i<FFT_LEN;i++){
            HX1[i].real = H[i]*X1[i].real;
            HX1[i].imag = H[i]*X1[i].imag;
        }
        //TODO ifft
        ifft(HX1, FFT_LEN);
        for(i=0;i<FFT_LEN;i++){
            enhSpeech_Frame_tmp[i] = HX1[i].real;
        }
        for(i=0;i<frameLength;i++){
            enhSpeech_Frame[i] = enhSpeech_Frame_tmp[i];
        }
        for(i=iniFrameSample-1;i<frameLength;i++){
            enhanced_output[i] = enhanced_output[i]+enhSpeech_Frame[i];
        }
        //update frame boundaries
        iniFrameSample = iniFrameSample + frameShift;
        endFrameSample = endFrameSample+frameShift;
    }

    //free array allocations
    //free window
    free(window);
    //free enhanced_output
    for(i=0;i<lenS;i++){
        ptrFloat = enhanced_output[i];
        free(ptrFloat);
    }
    //free P
    for(i=0;i<FFT_LEN/2;i++){
        ptrFloat = P[i];
        free(ptrFloat);
    }
    //free limNeg
    for(i=0;i<FFT_LEN/2;i++){
        ptrFloat = limNeg[i];
        free(ptrFloat);
    }
    //free Frame1 and Frame2
    free(Frame1);
    free(Frame2);
    //free wFrame1 and wFrame2
    free(wFrame1);
    free(wFrame2);
    //free X1 and X2
    free(X1);
    free(X2);
    //free PX1X1 PX2X2 PX1X2
    free(PX1X1);
    free(PX2X2);
    free(PX1X2);
    //free cohX
    free(cohX);
    //free reCOH and imCOH
    free(reCOH);
    free(imCOH);
    //free G1 G2
    for(i=0;i<FFT_LEN/2;i++){
        ptrFloat = G1[i];
        free(ptrFloat);
    }
    for(i=0;i<FFT_LEN/2;i++){
        ptrFloat = G2[i];
        free(ptrFloat);
    }
    //free ind_neg
    free(ind_neg);
    //free H
    for(i=0;i<FFT_LEN/2;i++){
        ptrFloat = H[i];
        free(ptrFloat);
    }
    //free temporary variable HX1
    free(HX1);
    //free enhSpeech_Frame_tmp and enhSpeech_Frame
    free(enhSpeech_Frame);
    free(enhSpeech_Frame_tmp);
}