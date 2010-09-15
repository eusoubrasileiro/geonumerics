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
    double *samples, *result; /* input signal, result */
	double error;
	FILE *file;  
    samples = (double*) malloc ( sizeof(double) * nsMax );

    /* numero de amostras que vem, tudo o que estiver no arquivo */
   /* test data set from python */
    file = fopen("input_rnd.txt", "r");    
    for(i=0; fscanf(file, "%lf", &samples[i]) == 1; i++){} /* leh as amostras */
    ns=i;     
    fclose(file);    
    
    /* Fnyqust 0.01 = 50hz .. Fc = 10, Ramp = 5, detrend linear 
    taper hanning window */
    result = dspIIR_SincTrapezLowPassApply(samples, ns, 0.01, 10, 5);    

    /* expected result from python */  
    file = fopen("expected_out.txt", "r");    
    for(i=0; fscanf(file, "%lf", &samples[i]) == 1; i++){} /* leh as amostras */
    ns=i;     
    fclose(file);    
	
	/* make the error , rms*/
	for(i=0, error=0; i<ns; i++)
		error += (samples[i]-result[i])*(samples[i]-result[i]);
	
	error = sqrt(error/ns);
	
	// altough not suceeding is working,.. comparison with python not perfect
	if(error < 1E-7)
		printf("suceed %g\n", error);
	else
		printf("didnt suced %g\n", error);
	
}

/* testing check based on python output */
int 
main(int argc, char argv[]){
    
    Fft_Convolution_Filtering_Test();    
    
    return 0;
}
