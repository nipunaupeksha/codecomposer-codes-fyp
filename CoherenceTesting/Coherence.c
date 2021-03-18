#include "stdio.h"
#include <math.h> /*Definitions for math library*/

#define PI 3.1415926535897

float* hanning_window(int frame_length);

/* fs - sampling frequency
 * x1 - signal 1
 * x2 - signal 2*/
void coherence(int fs, int* x1, int* x2 ){
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

	
	//float window[frameLength] = {0};
	//hanning_window(window,frameLength);
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
