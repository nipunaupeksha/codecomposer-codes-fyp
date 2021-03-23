#define M_PI        3.14159265358979323846      // Pi

void voiceActivityDetector(double* input, int fs){
    int i;//iterators
    double* inputSignal;double* vadMarker; double* outputSignal;
    int frameLength; int frameLengthSamples; int numFrames; int additionalSamples;
    double* additionalArray;
    int inputSignalLength;//additional variables
    double segmentedSignal;//TODO: currently declare this as a double. Ask from chathuki about it.
    double* iVAD; double* fVAD;
    int hangoverThreshold; int inactiveFrameCounter;int firstVoicedFlag;
    int applyFadeOut; int bufferFilling; int bufferLimit; int bufferFull;double* bufferVoicedFrames;
    int energyMin; int energyMax; int deltaEmin; int deltaEmax; int lambda; int threshold;
    int initialEnergyMin; int initialNoiseDuration; int lowerEnergyVoiceBand; int upperEnergyVoiceBand;
    int periodicityThreshold;int energyRatioThreshold;
    double* RMSEArray; double* periodicityArray; double* ratioArray; double* thresholdArray; double* thresholdMarker;
    int fadeInTime; int fadeOutTime;
    double* fadeInEnvelope; double* fadeOutEnvelope;
    double *frame; int frameSize; double RMSE;
    double periodicity;


    frameLength = 10; //in ms
    frameLengthSamples = ceil((double)(frameLength*fs)/1000.0); //no. of samples in one frame
    numFrames = ceil((double)inputSignalLength/frameLengthSamples);
    additionalSamples = (frameLengthSamples*numFrames)-inputSignalLength;
    additionalArray = zeros1d(additionalSamples); // no. of addtional samples in all frames
    inputSignal = (double*)calloc(inputSignalLength+additionalSamples, sizeof(double));
    //now input signal has more samples
    for(i=0;i<(inputSignalLength+additionalSamples);i++){
        if(i<inputSignalLength){
            inputSignal[i] = input[i];
        }else{
            inputSignal[i] = additionalArray[i];
        }
    }

    //Segment the signal into frames
    segmentedSignal = frameSplit(inputSignal, frameLengthSamples); //TODO: I defined frameSplit method in frameSplit.c  as a void function. But that must be changed.

    //Initialize VAD flag arrays
    iVAD = zeros1d(numFrames);
    fVAD = zeros1d(numFrames);

    //Hang over smoothing variables
    hangoverThreshold = 5;
    inactiveFrameCounter = hangoverThreshold +1;
    firstVoicedFlag = 0;
    applyFadeOut = 0;
    bufferFilling = 0;
    bufferLimit = 4;
    bufferFull = 0;

    //Initialize Energy variables
    energyMin = 0;
    energyMax = 0;
    deltaEmin = 1;
    deltaEmax = 1;
    lambda = 1;
    threshold = 0;
    initialEnergyMin = 0;
    initialNoiseDuration = 50;
    lowerEnergyVoiceBand = 2000;
    upperEnergyVoiceBand = 4000;

    //Initialize Periodicity and Energy Ratio Thresholds
    periodicityThreshold = 1;
    energyRatioThreshold = 10;

    //Initilize RMSE, Periodicity, and energy ratio arrays
    RMSEArray = zeros1d(numFrames);
    periodicityArray =zeros1d(numFrames);
    ratioArray = zeros1d(numFrames);
    thresholdArray = zeros1d(numFrames);

    //Set the fade in and fade out times for hangover smoothing
    fadeInTime = frameLengthSamples;
    fadeOutTime = frameLengthSamples*hangoverThreshold;

    //use a simple linear envelope for the fades
    fadeInEnvelope = linspace(0,1,fadeInTime);
    fadeOutEnvelope = linspace(0,1,fadeOutTime);

    //compute the thresholds and enfore the VAD process per frame
    for(i=0;i<numFrames+1;i++){
        //TODO: frame = segementedSignal(i,:) -->I guess its a 2D array
        //Get the Root Mean Square Energy of the frame
        RMSE = sqrt(mean(sum(frame), frameSize));

        //NOTE: An important assumption necessarily made is that the first 100ms
        //or so of the signal will be noise. So we use the first few frames to calculate
        // and set Emin and Emax values.
        if((frameSize*frameLength)<initialNoiseDuration){
            if(RMSE > energyMax){
                //set Emax
                energyMax = 5*RMSE;
            }
            if(energyMin == 0){
                //Initialize Emin
                energyMin = RMSE;
                initialEnergyMin = energyMin;
            }else{
                if(RMSE>energyMin){
                    energyMin = RMSE;
                }
            }
            iVAD[i] = 0;
            fVAD[i] = 0;
        }else{
            //Commence actual VAD processing
            //Updare running Energy threshold estimates(Emax, Emin)
            if(RMSE>energyMax){
                //We use a small delta to gradually decrease Emax to compensate
                //for anomalous spikes in enery
                energyMax = RMSE;
                deltaEmax = 1;
            }else{
                deltaEmax = 0.999;
            }

            if(RMSE<energyMin){
                if(RMSE==0){
                    energyMin = initialEnergyMin;
                }else{
                    energyMin = RMSE;
                }
                deltaEmin=1;
            }else{
                //Threshold computation. Lambda is the non-linear dynamic coefficient used to compute the threshold
                //in a way that makes it resistant and independent of variations in background noise.
                lambda = 1 - ((double)energyMin/energyMax);
                threshold = (((1-lambda)*energyMax) + (lambda*energyMax));

                //Keep track of threshold values
                thresholdArray[i] = threshold;

                //Get the periodicity of the frame
                periodicity = peakPeriodicty(frame,fs);

                //Add to periodicity array
                periodicityArray[i] = periodicity;

                //Get the ratio of the frequencies above and below 2kHz in the voice band(0-4kHz). We will then use the ratio
                //of these energies to make decision around the voicing of the frame.
                
                //TODO: window = hamming(size(frame, 2));

                //FFT computation and normalization
                //TODO: fftLength = 2^nextpow2(length(frame));
                //TODO: theFFT = fft(frame.*window, fftLength);
                //TODO: fftLength = length(theFFT);
                //TODO: fftSq = (theFFT).*conj(theFFT);
                //TODO: fftSqNorm = fftSq/fftLength;
            }
        }
    }

    //free callocs
    free(inputSignal);
    free(additionalArray);
    free(iVAD);
    free(fVAD);
    free(RMSEArray);
    free(bufferVoicedFrames);
    free(periodicityArray);
    free(ratioArray);
    free(thresholdArray);
    free(thresholdMarker);
    free(fadeInEnvelope);
    free(fadeOutEnvelope);
}

double* zeros1d(int n){//1D zeros vector
    double* z;
    z = (double*)calloc(n, sizeof(double));
    return z;
}

double* linspace(double start, double stop, int n){
    int i; //iterator
    double *vector;
    double step = (double)(stop - start)/n;
    vector = (double*)calloc(n,sizeof(double));
    for(i=0; i<n; i++){
        if(i!=n-1){
            vector[i] = start+step*i;
        }else{
            vector[i] = stop;
        }
    }
    return vector;
}

double sum(double* arr){
    double total;
    int arrLength = *(&arr+1)-arr;
    int i;
    for(i=0;i<arrLength;i++){
        total += arr[i];
    }
    return total;
}

double mean(double total, double size){
    return total/size;
}

double* Hamming(int N){
    double* w;
    int i;//Iterator
    w = (double*)calloc(N, sizeof(double));
    for(i=0;i<N;i++){
        w[i] = 0.54 - 0.46*cos((double)(2*M_PI*i)/(N-1));
    }
    return w;
}

double nextpow2(double n){
    return ceil(log2(abs(n)));
}