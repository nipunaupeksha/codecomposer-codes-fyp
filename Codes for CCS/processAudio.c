void processAudio(int* clean_speech, int* processed_speech){
    int i;
    int clear_speech_length; int processed_speech_length;
    double mc,mp; //mean
    double* in;
    double maxcs; double maxps;
    clear_speech_length = *(&clean_speech+1)-clean_speech;
    processed_speech_length = *(&processed_speech+1)-processed_speech;
    mc = mean(clean_speech,clear_speech_length);//mean of clear signal
    mp = mean(processed_speech,processed_speech_length);//mean of processed signal
    for(i=0;i<clear_speech_length;i++){
        clean_speech[i] = clean_speech[i] - mc; 
    }
    for(i=0;i<processed_speech_length;i++){
        processed_speech[i] = processed_speech_length - mp;
    }
    in = (double*)calloc(processed_speech_length, sizeof(double));
    maxcs = absMax(clean_speech, clear_speech_length);
    maxps = absMax(processed_speech,processed_speech_length);
    for(i=0;i<processed_speech_length;i++){
        in[i] = processed_speech[i] * ((double)maxcs/maxps);
    }
    //free callocs
    free(in);
}

double mean(double* arr, int length){
    int i;
    double total;
    for(i=0;i<length;i++){
        total+=arr[i];
    }
    return (double)total/length;
}

double absMax(double* arr, int length){
    int i;
    double maxValue;
    maxValue = 0;
    for(i=0;i<length;i++){
        maxValue = maxValue<abs(arr[i])? abs(arr[i]): maxValue; 
    }
    return maxValue;
}