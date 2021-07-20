#ifndef COHERENCE_H_
#define COHERENCE_H_

/*definitions*/
#define lambdaX 0.68 //forgetting factor for smoothing power spectrum
#define epsilon 10^(-12)
#define negImgConsFilter 0.05

float* COH_RealImag(Int16 *x1, Int16 *x2, Int32 fs);

#endif /*COHERENCE_H_*/
