void frameSplit(int *input, int windowLength)
{
    int i, j; //iterator
    float *finalSignal;
    float *indexFactor;
    int *indices;
    int numFrames;
    int inputSize;
    int increment;
    inputSize = *(&input + 1) - input;

    increment = windowLength;

    numFrames = fix((inputSize - windowLength + increment) / (float)increment);
    indexFactor = (float *)calloc(numFrames, sizeof(float));
    memset(indexFactor, 0, sizeof(float));
    for (i = 0; i < numFrames; i++)
    {
        indexFactor[i] = increment * i;
    }
    for (i = 0; i < windowLength; i++)
    {
        indices[i] = i;
    }

    finalSignal = zeros2d(numFrames, windowLength);
    for (i = 0; i < numFrames; i++)
    {
        for (j = 0; j <= windowLength; j++)
        {
            //TODO: finalSignal(:) = input(indexFactor(:,ones(1,len))+indices(ones(numFrames,1),:)); %numFrames*len matrix
            finalSignal[i][j] = i * windowLength + indices[j]; 
            //0                         0                   0              ...             0
            //windowLength          windowLength        windowLength       ...        windowLength
            //2*windowLength        2*windowLength     2*windowLength      ...      2*windowLength
            //...
            //(windowLength)^2    (windowLength)^2    (windowLength)^2     ...  (windowLength)^2

        }
    }

    free(indexFactor);
    free(finalSignal);
    /* Construct final signal -> since windowLength is not a matrix we dont need to implement this
    if (winLength > 1)
        winFinal = windowLength(:)';
        finalSignal = finalSignal .* winFinal(ones(numFrames,1),:);*/
}

int fix(float n)
{ //n =  number
    return n <= 0 ? (int)(ceil(n)) : (int)(floor(n));
}

float *zeros1d(int n)
{ //n = number
    float *z;
    z = (float *)calloc(n, sizeof(float));
    memset(z, 0, sizeof(float));
    return z;
}

float *zeros2d(int r, int c)
{ //r = rows, c = columns
    float *z;
    z = (float *)malloc(r * c * sizeof(float));
    int i, j;
    for (i = 0; i < r; i++)
    {
        for (j = 0; j < c; j++)
        {
            *(z + i * c + j) = 0;
        }
    }
    return z;
}

float *ones(int n)
{ //n = number
    float *o;
    o = (float *)calloc(n, sizeof(float));
    memset(o, 1, sizeof(float));
    return o;
}