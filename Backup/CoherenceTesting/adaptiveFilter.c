#include "conf.h"
#include <math.h>

void adaptive_filter(Int32* x){
	Int32 N; /*Length of noisy input signal*/
	
	N = *(&x+1)-x;
	
}