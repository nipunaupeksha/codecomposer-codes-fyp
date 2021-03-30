#include "stdio.h"
#include "conf.h"
#include <math.h>

#define lambda_x 0.68 //Forgetting factor for smoothing power specturm
#define epsilon 10^(-12)
#define negImag_ConsFilter 0.05 //Constant value of filter when imag part is negative
// 
//   x1, x2 : Input signals at two channels
//   fs : Sampling frequency
//   
void COH_RealImag(double* x1, double* x2, Int32 fs){
	Int32 i; //iterators
	Int32 frameLength, frameShift, FFT_LEN;
    float* window; float* enhanced_output; 
    int lenS, nFrame, iniFrameSample, endFrameSample;
	float band1;
	float* P; 
    float* limNeg; 
    float* Frame1; float* Frame2; 
	float* wFrame1; float* wFrame2; 
    float* PX1X1; float* PX2X2; complex* PX1X2;  
    float* reCOH; float* imCOH; 
	float* G1; float* G2; float* onesG; float* G; 
    float* H1; float* H2;
    complex* v1; complex* v2; complex* scratch; 
    complex* cohX; complex* enhSpeech_Frame; 

	
	frameLength = (Int32)floor(fs*20/(float)1000);
	if((frameLength/(float)2) != floor((frameLength/(float)2))){
		frameLength += 1;
	}
	
		
	frameShift = (Int32)floor(frameLength * 0.5);
	window = hanning(frameLength);
	FFT_LEN = (Int32)exp2(frameLength);
	
	lenS = min((*(&x1+1)-x1),(*(&x2+1)-x2));
	nFrame = 0;
	iniFrameSample = 0; //I think it should be 0. But in MATLAB code it is initialized as 1
	endFrameSample = iniFrameSample + frameLength -1;
	enhanced_output = zeros(lenS);
	
	//defining bands
	band1 = floor(1000*FFT_LEN/(float)fs); //0 -> 1 kHz
	//Defining exponents of coherence function
	
	P = zeros(FFT_LEN/2);
	//Defining thresholds for imaginary parts to consider negative(noise)
	
	limNeg = zeros(FFT_LEN/2);
	
	for(i=0;i<band1;i++){
		P[i] = 16;
		limNeg[i] = -0.1;
	}
	for(i=band1;i<FFT_LEN/2;i++){
		P[i] = 2;
		limNeg[i] = -0.3;
	}
	
	//set sizes of arrays for processing
	Frame1 = (float*)calloc(frameLength,sizeof(float));
	Frame2 = (float*)calloc(frameLength,sizeof(float));
	
	wFrame1 = (float*)calloc(frameLength,sizeof(float));
	wFrame2 = (float*)calloc(frameLength,sizeof(float));
	
	v1 = (complex*)calloc(frameLength, sizeof(complex));
	v2 = (complex*)calloc(frameLength, sizeof(complex));
    scratch = (complex*)calloc(frameLength, sizeof(complex));

    PX1X1 = (float*)calloc(frameLength, sizeof(float));
    PX2X2 = (float*)calloc(frameLength, sizeof(float));
    PX1X2 = (complex*)calloc(frameLength, sizeof(complex));

    G1 = (float*)calloc(FFT_LEN/2, sizeof(float));
    G = (float*)calloc(FFT_LEN/2, sizeof(float));
    enhanced_output = (float*)calloc(sizeof(float));

    H1 = (float*)calloc(FFT_LEN/2, sizeof(float));
    H2 = (float*)calloc(FFT_LEN/2, sizeof(float));

    reCOH = (float*)calloc(FFT_LEN/2,sizeof(float));
    imCOH = (float*)calloc(FFT_LEN/2,sizeof(float));
    cohX = (complex*)calloc(FFT_LEN/2,sizeof(float));
    enhSpeech_Frame = (complex*)calloc(frameLength, sizeof(complex));

	while(endFrameSample<lenS){
		nFrame += 1;//A new frame to process
		
		//Get short-time magnitude and phase spectrum for each input channel
		for(i=0;i<frameLength;i++){
    		Frame1[i] = x1[iniFrameSample+i];
    		Frame2[i] = x2[iniFrameSample+i];
		}
		for(i=0;i<frameLength; i++){
			wFrame1[i] = Frame1[i]*window[i];
			wFrame2[i] = Frame2[i]*window[i];
		}
		
		for(i=0;i<frameLength;i++){
			v1[i].real = wFrame1[i];
			v1[i].imag = 0;
			v2[i].real = wFrame1[i];
			v2[i].imag = 0;
			scratch[i].real =0;
			scratch[i].imag = 0;
		}
		fft(v1,FFT_LEN,scratch);
		fft(v2, FFT_LEN,scratch);
		if(nFrame == 1){
			for(i=0;i<frameLength;i++){
				PX1X1[i] = (v1[i].real * v1[i].real) + (v1[i].imag * v1[i].imag);
				PX2X2[i] = (v2[i].real * v2[i].real) + (v2[i].imag * v2[i].imag);
				PX1X2[i].real = (v1[i].real * v2[i].real) + (v1[i].imag * v2[i].imag);
				PX1X2[i].imag = (v1[i].real * v2[i].imag) + (v1[i].imag * v2[i].real);
			}
		}else{
			for(i=0;i<frameLength;i++){
				PX1X1[i] = lambda_x * PX1X1[i] + (1-lambda_x)*(v1[i].real * v1[i].real) + (v1[i].imag * v1[i].imag);
				PX2X2[i] = lambda_x * PX2X2[i] + (1-lambda_x)*(v2[i].real * v2[i].real) + (v2[i].imag * v2[i].imag);
				PX1X2[i].real = lambda_x * PX1X2[i].real + (1-lambda_x)*(v1[i].real * v2[i].real) + (v1[i].imag * v2[i].imag);
				PX1X2[i].imag = lambda_x * PX1X2[i].imag + (1-lambda_x)*(v1[i].real * v2[i].imag) + (v1[i].imag * v2[i].real);
			}
		}
		for(i =0;i<frameLength;i++){
			cohX[i].real = PX1X2[i].real / (float)(PX1X1[i]*PX2X2[i]);
			cohX[i].imag = PX1X2[i].imag / (float)(PX1X1[i] * PX2X2[i]);
		}
		for(i=2;i<=FFT_LEN/2+1;i++){
			reCOH[i-2] = cohX[i].real;
			imCOH[i-2] = cohX[i].imag;
		}
		onesG = ones(FFT_LEN/2);
		for(i=0;i<FFT_LEN/2;i++){
			G1[i] = onesG[i] - (float)abs(reCOH[i])*P[i]; /*for suppressing noise coming from angles greater than 90*/
		}
		G2 = ones(FFT_LEN/2);
		for(i=0;i<FFT_LEN/2;i++){
			if(imCOH[i]<limNeg[i]){
				G2[i] = negImag_ConsFilter;
			}
		}
		/*Halfband final filter*/
		for(i=0;i<FFT_LEN/2;i++){
			G[i] = G1[i] * G2[i];
		}
		/*Fullband final filter*/
		H2 = flipud(G);
		for(i=0;i<FFT_LEN/2;i++){
			H1[i] = (float)abs(G[i]);
			H2[i] = (float)abs(H2[i]);
		}
		for(i=0; i<FFT_LEN; i++){
			if(i<FFT_LEN/2){
				v1[i].real = v1[i].real * H1[i];
				v1[i].imag = v1[i].imag * H1[i];
			}else{
				v1[i].real = v1[i].real * H2[i];
				v1[i].imag = v1[i].imag * H2[i];
			}
		}
		ifft(v1, FFT_LEN, scratch);
		for(i=0; i<frameLength; i++){
			enhSpeech_Frame[i].real = v1[i].real; //real is enough??
		}
		for(i = 0;i<frameLength;i++){
			enhanced_output[iniFrameSample+i] = enhSpeech_Frame[i].real + enhanced_output[iniFrameSample+i];
		}
		/*Update frame boundaries*/
		iniFrameSample+=frameShift;
		endFrameSample += frameShift;
	}
	/*Free variables*/
	free(window);
	free(limNeg);
	free(P);
	free(Frame1);
	free(Frame2);
	free(wFrame1);
	free(wFrame2);
	free(PX1X1);
	free(PX2X2);
	free(reCOH);
	free(imCOH);
	free(G1);
	free(G2);
	free(G);
	free(H1);
	free(H2);
	free(enhSpeech_Frame);
	free(PX1X2);
	free(v1);
	free(v2);
	free(scratch);
	free(cohX);
	
	return;
}	