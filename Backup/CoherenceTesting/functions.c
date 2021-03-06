#include "conf.h" /*Definition of complex variable structure*/
#include "stdlib.h"
#include "string.h"
#include <math.h> /*Definitions of math library*/

#define PI 3.14159265358979323846264338327950288
#define clipThresholdFactor 0.5
#define lowerLagMS 2.5
#define upperLagMS 5

/************************************************
 *   Implementation of Hann function in MATLAB
 * *********************************************/
float* hanning(Int32 N){
	Int32 i;
	float* w = (float*)malloc(sizeof(float) * N);
	memset(w, 0, N*sizeof(float));
	for(i = 0; i<N; i++){
		w[i] = 0.5 * (1 - cos(2 * PI * i/(float)(N-1)));
	}
	return w;
}

/************************************************
*    Implementation of Zeors function in MATLAB
*    1D Vector
************************************************/
float* zeros(Int32 N){
	float* z = (float*) malloc(N * sizeof(float));
	memset(z, 0, N*sizeof(float));
	return z;
}

/************************************************
*    Implementation of Ones function in MATLAB
*    1D Vector
************************************************/
float* ones(Int32 N){
	int i; //iterator
	float* o = (float*) malloc(N * sizeof(float));
	for(i-0;i<N;i++){
		o[i] = 1;
	}
	return o;
}

/*************************************************
*    Implementatiion of nextpow2 function in
*    MATLAB
*************************************************/

Int32 nextpow2(Int32 n){
	Int16 nextPow;
	n = abs(n);
	nextPow = ceil(log2(n));	
	return nextPow;
}

/*************************************************
 *    Implementation of flipup function in 
 *    MATLAB
 ************************************************/
 float* flipud(float* arr){
 	/*arr->an array*/
 	Int32 i; //Iteration
 	/* length of array
 	 * *(&arr + 1) - arr;
 	 * sizeof(arr)/sizeof(arr[0]) 
 	 */
 	Int32 N = *(&arr + 1) - arr;
 	float *f =(float*) malloc(N * sizeof(float));
 	memset(f, 0, N*sizeof(float));
 	for(i=N-1;i>=0;i--){
 		f[N-1-i] = arr[i];
 	}
 	return f;
 }

/*************************************************
*    Implementation of min function
*************************************************/
double min(double n1, double n2){
	return n1 <= n2 ? n1 : n2;
}

/************************************************
*    Implementation of max function
************************************************/
double max(double n1,double n2){
	return n1 >= n2 ? n1 : n2;
}

/*************************************************
 *   Periodicity Function -> Peak Periodicity
 ************************************************/
 void peakPeriodicity(Int32* input, Int32 fs){
 	Int32 i; //Iterator
 	Int32 duration;//length of input signal
 	float clipThreshold; float lowerLag; float upperLag;
 	float* output; float* acf;float* hVector;float* lagArray;
 	float maxValue=0;
 	
 	duration = *(&input)-input;
 	
 	for(i=0;i<duration;i++){
 		maxValue = (float)max(maxValue,input[i]);
 	}
 	
 	clipThreshold = clipThresholdFactor * maxValue;
 	output = zeros(duration);
	
	for(i=0;i<duration;i++){
		if(input[i]>=clipThreshold){
			output[i] = input[i] - clipThreshold;
		}else if(abs(input[i])<clipThreshold){
			output[i]= 0;
		}else if(input[i]<-clipThreshold){
			output[i] = input[i]+clipThreshold;
		}
	} 	
 	
 	//TODO acf = xcorr(output)
 	lowerLag = (float)round((lowerLagMS * fs)/(float)1000);
 	upperLag = (float)round((upperLagMS * fs)/(float)1000);
 	//TODO acfSectionalized = acf(duration+lowerLag:duration+upperLag)
 	for(i=0;i<upperLag-lowerLag+1;i++){
 		hVector[i] = lowerLag+i;
 	}
 	for(i=0;i<upperLag-lowerLag+1;i++){
 		lagArray[i] = 1 - hVector[i];
 	}
 	//TODO acfNormalized = acfSectionalized ./ lagArray;

	//TODO energyNormalized = acf(duration)/duration;

	//TODO finalacf = acfNormalized / energyNormalized;

	//TODO peakPeriodicity = max(finalacf);
 	
 }

/*************************************************
*    Voice Activity Detector
*************************************************/
/*
void voice_activity_detector(Int32* x, Int32 fs){
	Int32 i; //Iterator
	float* vadMarker; float* outputSignal; float* additionalArray;
	float* inputSignal;
	Int32 frameLength; Int32 frameLengthSamples; Int32 numFrames; Int32 signalLength;
	Int32 additionalSamples;
	
	
	signalLength = *(&x+1)-x;
	
	frameLength = 10;
	frameLengthSamples = (Int32)ceil(frameLength*fs/(float)1000);
	numFrames = (Int32)ceil(signalLength/(float)frameLengthSamples);
	additionalSamples = (frameLengthSamples*numFrames - signalLength);
	additionalArray = zeros(additionalSamples);
	for(i=0;i<signalLength;i++){
		
	} //TODO line 16/17
	
}*/

