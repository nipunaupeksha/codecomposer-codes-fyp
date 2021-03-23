void peakPeriodicity(int* input, int fs){
    int i; //iterators
    const float clipThreshodFactor = 0.5;
    float clipThreshold;
    float maxValue;
    int duration;
    float* output; //input and output are floats or ints
    /* TODO: I think I should use malloc function for all the pointer declarations. 
    /Lets see that after implementing this code*/
    double* acf; float lowerLagMS; float upperLagMS; int lowerLag; int upperLag;
    double* acfSectionalized; double* hVector; double* lagArray; double* acfNormalized;
    double* finalacf; double energyNormalized; double peakPeriodicity;
    maxValue = 0;
    duration = *(&input+1) - input;
    for(i=0;i<duration;i++){
        maxValue = max(input[i], maxValue);
    }
    clipThreshold = clipThreshodFactor * maxValue;

    //center clipping algorithm
    for(i=0;i<duration;i++){
        if(input[i] >=clipThreshold){
            output[i] = input[i] -clipThreshold;
        }else if(abs(input[i])<clipThreshold){
            output[i] = 0;
        }else if(input[i] < -clipThreshold){
            output[i] = input[i] + clipThreshold;
        }
    }

    acf = (double*)calloc(2*duration + 1,sizeof(double));
    xcorr(input, input, acf, duration, 0, 0, 0, 0, duration, duration, 0, 0);

    /*Now, we want the ACF only for lag values that fall within the pitch period limits. 
    * Let's take 2.5ms to 15ms, i.e. with a sampling rate of fs, that corresponds to*/
    lowerLagMS = 2.5;
    upperLagMS = 5;
    lowerLag = (float)round((lowerLagMS*fs)/(float)1000);
    upperLag = (float)round((upperLagMS*fs)/(float)1000);
    
    /*The ACF is symmetric around 0, meaning it goes from -n to n(lag values).
    * We want the ACF values for lags going from the lower to the upper lag limits.
    * Because there are 2n-1 ACF values, this maps to n+ the lag we want. */
    //TODO: acfSectionalized = acf(duration+lowerLag:duration+upperLag);
    acfSectionalized = (double*)calloc(upperLag - lowerLag + 1, sizeof(double));
    for(i=0;i<upperLag-lowerLag+1;i++){
        acfSectionalized[i] = acf[duration +lowerLag+i];
    }

    //We normalizing the ACF with n-h and then energy per sample
    hVector = (double*)calloc(upperLag-lowerLag+1, sizeof(double));
    lagArray = (double*)calloc(upperLag -lowerLag +1, sizeof(double));
    for(i=0;i<upperLag-lowerLag+1;i++){
        hVector[i] = lowerLag+i;
    }
    for(i=0;i<upperLag-lowerLag+1;i++){
        lagArray[i] = 1 - hVector[i];
    }
    acfNormalized = (double*)calloc(upperLag-lowerLag+1, sizeof(double));
    for(i=0;i<upperLag-lowerLag+1;i++){
        acfNormalized[i] = acfSectionalized[i]/(double)lagArray[i];
    }

    finalacf = (double*)calloc(upperLag-lowerLag+1, sizeof(double));
    energyNormalized = acf[duration-1]/duration;
    for(i=0;i<upperLag-lowerLag+1;i++){
        finalacf[i] = acfNormalized[i]/(double)energyNormalized;
    }
    peakPeriodicity = 0;
    for(i=0;i<upperLag-lowerLag+1;i++){
        peakPeriodicity = max(peakPeriodicity,finalacf[i]);
    }

    free(acf);
    free(acfSectionalized);
    free(hVector);
    free(lagArray);
    free(acfNormalized);
    //return peakPeriodicity if necessary
}

float max(float n1, float n2){
    return n1>=n2 ? n1 : n2;
}

float min(float n1, float n2){
    return n1<=n2 ? n1 : n2;
}

float* zeros1d(int n){ //n = number
    float* z;
    z = (float*)malloc(n * sizeof(float));
    memset(z,0, sizeof(float));
    return z;
}

void xcorr(float *tr1, float *tr2, double *corp, int shift, int shift_zero,
            int window, int dmean, int normalize, int ndat1, int ndat2, 
            int ndat1d, int ndat2d){
    // Calculates the cross-correlation function of data arrays tr1 and tr2.
    // We use the following definition for cross-correlation:
    // corp[i] = sum_j(tr1[i+j]*tr2[j])
    // The data is demeaned before. The result is normalized.
    // input:
    // tr1, tr2     : data arrays with length ndat1, ndat2
    // shift        : maximal shift
    // shift_zero   : before cross-correlation the first data array is right-shifted
    //                by that amount of samples.
    //                (The effect is a right-shift of the correlation function)
    //                In the formula above tr1[i+j] is substituted by tr1[i+j-shift_zero]
    // dmean        : 0 or 1. if 1 the data is demeand
    // normalize    : 0 or 1. if 1 the cross-correlation function is normalized
    //                (correlation coefficient of 1 than means perfect fit)
    // window       : only use values in this window for demeaning and normalization
    //                if 0  : window = min(ndat1, ndat2)
    //                if >0 : window = this parameter
    // ndat1, ndat2 : length of data arrays
    // ndat1d, ndat2d: use this length for demeaning(see source code)
    //                 (if 0 ndat1d = ndat2d = window)
    // output       : 
    // corp         : cross-correlation function of length 2*shift+1
    int a,a2,b,b2,bmin,bmax,flag=0, ind1,ind2,ind3,ind4;
    double sum, sum1, sum2, cmax;
    float *tra1, *tra2;
    tra1 = (float*)calloc(ndat1, sizeof(float));
    tra2 = (float*)calloc(ndat2, sizeof(float));

    //Set standard values
    if(window == 0){
        window = (int)min((float)ndat1, (float)ndat2);
    }
    if(ndat1d == 0){
        ndat1d = window;
    }
    if(ndat2d == 0){
        ndat2d = window;
    }
    ind1 = (int)max(0,(float)(ndat1-window)/2);
    ind2 = (int)min(ndat1,(float)(ndat1+window)/2);
    ind3 = (int)max(0,(float)(ndat2-window));
    ind4 = (int)min(ndat2,(float)(ndat2+window)/2);

    // Demean data (Zero offset)
    if(dmean>0){
        sum = 0;
        for(a = ind1; a<ind2; a++){
            sum += tr1[a];
        }
        sum /= ndat1d;
        for (a = 0; a < ndat1; a++) {
            tra1[a] = tr1[a] - (float) sum;
        }
        if (sum == 0.0){
            flag = 1;
        }
        sum = 0;
        for (a = ind3; a < ind4; a++) {
            sum += tr2[a];
        }
        sum /= ndat2d;
        for (a = 0; a < ndat2; a++) {
            tra2[a] = tr2[a] - (float) sum;
        }
        if (sum == 0.0){
            flag += 1;
        }
    }else{
        for (a = 0; a < ndat1; a++) {
            tra1[a] = tr1[a];
        }
        for (a = 0; a < ndat2; a++) {
            tra2[a] = tr2[a];
        } 
    }
     if (flag == 0) {
        a = 0;
        a2 = -shift_zero - shift;
        if (ndat1 != ndat2) {
            for (; a < (2 * shift + 1); a++, a2++) {
                bmin = max(0, -a2 + (ndat2 - ndat1) / 2);
                bmax = min(ndat2, -a2 + (ndat1 + ndat2) / 2);
                b2 = bmin;
                b = b2 + (ndat1 - ndat2) / 2;
                corp[a] = 0;
                //printf("%d - %d %d - %d %d\n", a2, b + a2, b + a2 + bmax - bmin, b2, bmax);
                if (bmin >= bmax) {
                    continue;
                }
                for (; b2 < bmax; b++, b2++) {
                    corp[a] += tra1[b + a2] * tra2[b2];
                }
            }
        } else { // same as above only ndat2 = ndat1, we need one variable less in the second loop
            for (; a < (2 * shift + 1); a++, a2++) {
                bmin = max(0, -a2);
                bmax = min(ndat1, -a2 + ndat1);
                corp[a] = 0;
                if (bmin >= bmax) {
                    continue;
                }
                for (b = bmin; b < bmax; b++) {
                    corp[a] += tra1[b + a2] * tra2[b];
                }
            }
        }

        /* normalize xcorr function */
        if (normalize > 0) {
            sum1 = sum2 = 0.0;
            for (a = ind1; a < ind2; a++) {
                sum1 += (*(tra1 + a))*(*(tra1 + a));
            }
            for (a = ind3; a < ind4; a++) {
                sum2 += (*(tra2 + a))*(*(tra2 + a));
            }
            sum1 = sqrt(sum1);
            sum2 = sqrt(sum2);
            cmax = 1 / (sum1 * sum2);
            for (a = 0; a < (2 * shift + 1); a++) {
                corp[a] *= cmax;
            }
        }
    } else {
        for (a = 0; a < (2 * shift + 1); a++) {
            corp[a] = 0;
            //shift_max = 0;
            //corp_max = 0;
        }
    }
    free((char *) tra1);
    free((char *) tra2);

}

/*
float xcorr(int *x, int *y){
    int i, j, delay; //itertors
    int n = *(&x+1)-x;
    float mx, my, sx, sy, sxy, denom, r;
    int delay_max = n;
    
    // x , y = input signals --> same length
    // mx = mean of x signal
    // my = mean of y signal
    // sx
    // sy
    // sxy
    // denom = denominator
    // r
    
   //calculate the mean of two series
   mx = 0;
   my = 0;
   for(i=0;i<n;i++){
       mx += x[i];
       my += y[i];
   }
   mx /= n;
   my /= n;

   //calculate the denominator
   sx = 0;
   sy = 0;
   for(i=0;i<n;i++){
       sx += (x[i]-mx)*(x[i]-mx);
       sy += (x[i]-my)*(x[i]-my);
   }
   denom = sqrt(sx * sy);

   //calculate the correlation series
   for(delay=-delay_max;delay<delay_max;delay++){
       sxy = 0;
       for(i=0;i<n;i++){
           j = i + delay;
           while(j<0){
               j+=n;
           }
           j%=n;
           sxy = (x[i] -mx)*(y[j]-my);
       }
       r = sxy/denom;
   }
    return r;
}
*/