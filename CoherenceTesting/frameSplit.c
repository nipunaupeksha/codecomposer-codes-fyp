#include <math.h>
#include "stdio.h"

void frameSplit(int* input, int windowLength, int increment){
	int inputSize; int winLength; int len; int numFrames;
	
	inputSize = *(&input+1)-input;
	winLength = sizeof(windowLength)/sizeof(int);
	
	if(winLength==1){
		len = windowLength;
	}else{
		len = winLength;
	}
	if(increment<=0){
		increment = len;
	}
	numFrames = fix((inputSize-len+increment)/(float)increment);  
}

int fix(int n){
	if(n>=0){
		return floor(n);
	}else{
		return ceil(n);
	}
}