#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "..\c\dsprocessing.c"

/*
Main test for Fft_Convolution_Filtering_Demo, files or whatever as input (stdin)
whatch out for the maximum size
*/
#define nsMax 100000
void
Fft_Convolution_Filtering_Test(){
    unsigned int i=0, ns; /* counter, input size */
    double* samples; /* input signal */
	FILE *file;  
    samples = (double*) malloc ( sizeof(double) * nsMax );

    /* numero de amostras que vem, tudo o que estiver no arquivo */
  
    file = fopen("rnd.txt", "r");
    
    for(i=0; fscanf(file, "%lf", &samples[i]) == 1; i++){} /* leh as amostras */
    ns=i;  
    
    fclose(file);
    
	vPrintf(stdout, samples, ns);
    samples = dspIIR_SincTrapezLowPassApply(samples, ns, 0.01, 10, 5);    
    vPrintf(stdout, samples, ns);
}

int 
main(int argc, char argv[]){
    Fft_Convolution_Filtering_Test();    
    
    return 0;
}
