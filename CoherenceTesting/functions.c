#include "conf.h" /*Definition of complex variable structure*/
#include <math.h> /*Definitions of math library*/

#define PI 3.14159265358979323846264338327950288


/************************************************
 *   Implementation of Hann function in MATLAB
 * *********************************************/
float* hann(Int32 N){
	Int32 i;
	float* w;
	w = (float*) calloc(N, sizeof(float));
	for(i = 0; i<N; i++){
		w[i] = 0.5 * (1 - cos(2 * PI * i/(double)(N-1)));
	}
	return w;
}

/************************************************
*    Implementation of Zeors function in MATLAB
*    1D Vector
************************************************/
float* zeros(Int32 N){
	float* z;
	z = (float*)calloc(N, sizeof(float));
	return z;
}