#include "stdio.h"
#include "conf.h"
#include <math.h> /*Definitions for math library*/
//#include <complex.h> /*Definitions for complex library*/

#define PI 3.1415926535897


float* hanning_window(int frame_length);
float* zeros(int length);
//void _fft(cplx buf[], cplx out[], int n, int step);
//void fft(cplx buf[], int n);

/* fs - sampling frequency
 * x1 - signal 1
 * x2 - signal 2*/
void coherence(int fs, int* x1, int* x2 ){
	const float negImag_ConsFilter = 0.05;
	int frameLength;
	float lambda_x=0.68;   /*Forgetting factor for smoothing power spectrum*/ 
	float epsilon=10^-12;
	float *window;
	float shift;
	int lenS;
	int nFrame =0;
	int iniFrameSample = 1;
	int endFrameSample;
	int FFT_LEN;
	float *enhanced_output;
	float *P;
	float band1;
	float *limNeg;
	int *Frame1;
	int *Frame2;
	float *wFrame1;
	float *wFrame2;
	float *X1;
	float *X2;
	float *PX1X1;
	float *PX2X2;
	complex *PX1X2;
	complex	*cohX;
	
	/*Iteration Values*/
	int i;
	
	fs = 16000; /*16 kHz sampling rate*/  /*TODO-->should be coming first I guess*/
	frameLength = floor(20 * fs/ 1000); /*Frame length*/
	
	/*If integer division of frameLength/2 is not equal to float division with floor function then add 1 to frameLength*/
	if(frameLength/2 != floor(frameLength/2)){ 
		frameLength = frameLength + 1;
	}
	
	shift = floor(frameLength * 0.25);
	window = hanning_window(frameLength);
	
	FFT_LEN = (int)exp2(ceil(log2(frameLength))); /*MATLAB nextpow2 -> here log2 gets the power of the value frameLength, then ceil function gets the next integer value*/
	
	lenS = sizeof(x1)>=sizeof(x2)?sizeof(x2):sizeof(x1);/*minimum length of signal x1 and x2*/	
	endFrameSample = iniFrameSample+frameLength-1;

	enhanced_output = zeros(lenS);
	band1 = floor(1000*FFT_LEN/fs);
	/*Defining exponents of coherence function*/
	P = zeros(FFT_LEN/2);
	for(i=0;i<band1;i++){
		P[i] = 16;
	}
	for(i=band1;i<FFT_LEN/2;i++){
		P[i] = 2;
	}
	
	/*Defining thresholds for imaginary parts to consider negative(noise)*/
	limNeg = zeros(FFT_LEN/2);
	for(i=0;i<band1;i++){
		limNeg[i] = -0.1;
	}
	for(i=band1;i<FFT_LEN/2;i++){
		limNeg[i] = -0.3;
	}

	while(endFrameSample<lenS){
		/*A new frame to process*/
		nFrame = nFrame+1;
		/*Get short time magnitude and phase spectrums for each input element*/
		for(i=iniFrameSample;i<endFrameSample+1;i++){
			Frame1[i] = x1[i];
		}
		for(i=iniFrameSample;i<endFrameSample+1;i++){
			Frame2[i] = x2[i];
		}
		for(i=0;i<sizeof(Frame1);i++){
			wFrame1[i] = Frame1[i]*window[i];
		}
		for(i=0;i<sizeof(Frame2);i++){
			wFrame2[i] = Frame2[i]*window[i];
		}
		/*TODO fft*/
		if(nFrame==1){
			for(i=0;i<sizeof(X1);i++){
				PX1X1[i] = exp(abs(X1[i]),2);
			}
			for(i=0;i<sizeof(X2);i++){
				PX2X2[i] = exp(abs(X2[i],2));
			}
				/*TODO conjugate*/
		}else{
			for(i=0;i<sizeof(X1);i++){
				PX1X1[i] = lambda_x*PX1X1[i] + (1-lambda_x)*exp(abs(X1[i]),2);
			}
			for(i=0;i<sizeof(X2);i++){
				PX2X2[i] = lambda_x*PX2X2[i] + (1-lambda_x)*exp(abs(X2[i]),2);
			}
			for(i=0;i<sizeof(X1);i++){
				PX1x2[i] = lambda_x*PX1X2[i]+(1-lambda_x)*X[i]; /*TODO-conjugate*/
			}
		}
		for(i=0;i<sizeof(PX1X2);i++){
			cohX[i] = PX1X2[i]/sqrt(PX1X1[i]*PX2X2[i]+epsilon);
		}
		/*get real and imaginary coherence values*/
	}
}

float* hanning_window(int frame_length){
	float* window = (float*)malloc(sizeof(float) * frame_length);
	int i=0;
	for(i=0; i<frame_length; i++){
		float term = i/(float)(frame_length-1);
		window[i] = 0.5 * (1 - cos((2*PI)*term));
	}
	return window;
}

float* zeros(int length){
	float* array = (float*)malloc(sizeof(float)*length);
	int i=0;
	for(i=0;i<length;i++){
		array[i] = 0;
	}
	return array;
}

