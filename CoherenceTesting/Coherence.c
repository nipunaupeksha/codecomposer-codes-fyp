#include "stdio.h"
#include "conf.h"
#include <math.h>

#define lambda_x 0.68
#define epsilon 10^(-12)
#define newImag_ConsFilter 0.05

void COH_RealImag(Int32* x1, Int32* x2, Int32 fs){
	Int32 i; /*Iteration values*/
	Int32 frameLength, frameShift, FFT_LEN, lenS, nFrame, iniFrameSample, endFrameSample;
	float* window; float* enhanced_output; float* P; float* limNeg;
	float band1;
	Int32 *Frame1; Int32 *Frame2;
	
	frameLength = (Int32)floor(fs*15/(float)1000);
	if((frameLength/(float)2) != floor((frameLength/(float)2))){
		frameLength += 1;
	}
	
	Frame1 = (Int32*)calloc(frameLength, sizeof(Int32));
	Frame2 = (Int32*)calloc(frameLength, sizeof(Int32));
		
	frameShift = (Int32)floor(frameLength * 0.5);
	window = hanning(frameLength);
	FFT_LEN = (Int32)exp2(frameLength);
	lenS = min((*(&x1+1)-x1),(*(&x2+1)-x2));
	nFrame = 0;
	iniFrameSample = 1; //or 0?
	endFrameSample = iniFrameSample + frameLength -1;
	enhanced_output = zeros(lenS);
	band1 = floor(1000*FFT_LEN/(float)fs);
	P = zeros(FFT_LEN/2);
	limNeg = zeros(FFT_LEN/2);
	for(i=0;i<band1;i++){
		P[i] = 16;
		limNeg[i] = -0.1;
	}
	for(i=band1;i<FFT_LEN/2;i++){
		P[i] = 2;
		limNeg[i] = -0.3;
	}
	
	while(endFrameSample<lenS){
		nFrame += 1;
		
	}
}	