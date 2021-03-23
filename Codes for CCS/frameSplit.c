void frameSplit(int* input, int windowLength, int increment){
    int i,j; //iterator
    float* finalSignal;
    float* indexFactor;
    int* indices;
    int numFrames; int inputSize;
    inputSize = *(&input+1)-input;
    if(increment<=0){
        increment = windowLength;
    }
    numFrames = fix((inputSize-windowLength+increment)/(float)increment);
    indexFactor = (float*)malloc(numFrames * sizeof(float));
    memset(indexFactor, 0 ,sizeof(float));
    for(i=0;i<numFrames;i++){
        indexFactor[i] = increment*i;
    }
    for(i=0;i<windowLength;i++){
        indices[i] = i;
    }

    free(indexFactor);
    finalSignal = zeros2d(numFrames,windowLength);
    for(i=0;i<numFrames;i++){
        for(j=0;j<windowLength;j++){
            //TODO: finalSignal(:) = input(indexFactor(:,ones(1,len))+indices(ones(numFrames,1),:)); %numFrames*len matrix

        }
    }

    
    /* Construct final signal -> since windowLength is not a matrix we dont need to implement this
    if (winLength > 1)
        winFinal = windowLength(:)';
        finalSignal = finalSignal .* winFinal(ones(numFrames,1),:);*/
}

int fix(float n){ //n =  number
    return n<=0 ? (int)(ceil(n)): (int)(floor(n));
}

float* zeros1d(int n){ //n = number
    float* z;
    z = (float*)malloc(n * sizeof(float));
    memset(z,0,sizeof(float));
    return z;
}

float* zeros2d(int r,int c){ //r = rows, c = columns
    float* z;
    z = (float*)malloc(r * c * sizeof(float));
    int i, j;
    for (i = 0; i <  r; i++){
        for (j = 0; j < c; j++){
             *(z + i*c + j) = 0;
        }
    }
    return z;
}

float* ones(int n){ //n = number
    float* o;
    o = (float*)malloc(n * sizeof(float));
    memset(o, 1, sizeof(float));
    return o;
}