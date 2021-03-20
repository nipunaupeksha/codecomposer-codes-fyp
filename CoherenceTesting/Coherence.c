#include "stdio.h"
#include "conf.h"
#include <math.h>

#define lambda_x 0.68
#define epsilon 10^(-12)
#define negImag_ConsFilter 0.05

void COH_RealImag(Int32* x1, Int32* x2, Int32 fs){
	Int32 i; /*Iteration values*/
	Int32 frameLength, frameShift, FFT_LEN, lenS, nFrame, iniFrameSample, endFrameSample;
	float* window; float* enhanced_output; float* P; float* limNeg;
	float band1;
	float* Frame1; float* Frame2; float *wFrame1; float *wFrame2; float* PX1X1; float* PX2X2; 
	complex* PX1X2;
	complex *v1; complex *v2; complex *scratch;
	complex * cohX;
	float * reCOH; float* imCOH; float *G1; float *G2; float *G; float* H1;
	
	frameLength = (Int32)floor(fs*15/(float)1000);
	if((frameLength/(float)2) != floor((frameLength/(float)2))){
		frameLength += 1;
	}
	
	wFrame1 = (float*)calloc(frameLength, sizeof(float));
	wFrame2 = (float*)calloc(frameLength, sizeof(float));
		
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
				PX1X1[i] = (v1[i].real * v1[i].real) + (v1[i].imag * v[i].imag);
				PX2X2[i] = (v2[i].real * v2[i].real) + (v2[i].imag * v2[i].imag);
				PX1X2[i].real = (v1[i].real * v2[i].real) + (v1[i].imag * v2[i].imag);
				PX1X2[i].imag = (v1[i].real * v2[i].imag) + (v1[i].imag * v2[i].real);
			}
		}else{
			for(i=0;i<frameLength;i++){
				PX1X1[i] = lambda_x * PX1X1[i] + (1-lambda_x)*(v1[i].real * v1[i].real) + (v1[i].imag * v[i].imag);
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
		for(i=0;i<FFT_LEN/2;i++){
			G1 = 1 - (float)abs(reCOH[i])*P[i];
		}
		G2 = ones(FFT_LEN/2);
		for(i=0;i<FFT_LEN/2;i++){
			if(imCOH[i]<limNeg[i]){
				G2[i] = negImag_ConsFilter;
			}
		}
		for(i=0;i<FFT_LEN/2;i++){
			G[i] = G1[i] * G2[i];
		}
		H1 = flipud(G);
		
		iniFrameSample+=frameShift;
		endFrameSample += frameShift;
	}
}	