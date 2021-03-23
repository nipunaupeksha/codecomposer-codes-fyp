//Delay compensate for the FIR filter
void delay_compensate(double* sig, double* filtered_sig, int fs, int N, int delay){ //is it int delay or double delay?
    int i;//iterator
    double* tn; double*tt; double* sn; double *sf;
    tn = (double*)calloc(N, sizeof(double));
    for(i=0;i<N;i++){
        tn[i] = (double)i/fs;
    }
    tt = (double*)calloc(N-delay, sizeof(double));
    for(i=0;i<N-delay+1;i++){
        tt[i] = tn[i];
    }
    sn = (double*)calloc(N-delay, sizeof(double));
    for(i=0;i<N-delay+1;i++){
        sn[i] = sig[i];
    }
    //release callocs
    free(tn);
    free(tt);
    free(sn);
}