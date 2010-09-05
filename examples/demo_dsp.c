#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "..\source\c\dsprocessing.c"

/*
Main example for FFT and DFT, files or whatever as input (stdin)
whatch out for the maximum size
*/
#define nsMax 10000
void
Dft_Fft_Demo(){
    unsigned int i=0, ns, ns_2; /* counter, input size, power 2 size */
    double samples[nsMax]; /* input signal */
    double *an, *bn; /* real and imaginary coeficients */
	double *samples_2;	/* padded with zeros signal */
	unsigned int time_;

    /* numero de amostras que vem, tudo o que estiver no arquivo */

	printf("Type the input signal just the values, after type any non number character\n");
    for(i=0; scanf("%lf", &samples[i])==1; i++){} /* leh as amostras */
    ns=i;  

	/* Comparing FFT with DFT, 
	so the input signal must be power2 size 
	because my function doesnt work for not power 2 size */
	if(isPower2(ns)==-1)
	{
		/* if its not get the next, realloc the array padding with zeros*/
		ns_2 = nextPower2(ns);
		
		samples_2 = (double*) malloc(sizeof(double)*ns_2);
		
		vCopiav(samples, samples_2, ns);

		/* padd the rest with zeros */
		for(i=ns; i<ns_2; i++)
			samples_2[i] = 0;

		/* ready */
	}
	/* after padding with zeros */
	fprintf(stderr,"\n n %u  ns_2 %u\n", ns, ns_2);
    vPrintf(stderr, samples_2, ns_2);

	/* Dft first */     
    dspDft_(samples_2, ns_2, &an, &bn); /* calcula os coeficientes real imaginarios*/
    fprintf(stderr, "\n complex coefs Dft\n");
    for(i=0; i<ns_2; i++) /* imprime as belezinhas */
        printf("%lf %lf\n", an[i], bn[i]);
    
    free(an);
    free(bn);		
    
	/* Fft after */
    dspFft_(samples_2, ns_2, &an, &bn);        
    fprintf(stderr, "\n complex coefs Fft\n");
    for(i=0; i<ns_2; i++) /* imprime as belezinhas */
        printf("%lf %lf\n", an[i], bn[i]);

   free(an);
   free(bn);
}


/*
Convolution games
Demo
Not efficient algorithm though working
*/

void
Convolution_Demo(){
    double *h, *x, *y; /* impulsiva, entrada, saida */
    int Nh, Nx; /* tamanho dos sinais */
    int t; /* counters */

    printf("Type the impulse response for the sistem, size and data\n");
    scanf("%d", &Nh);
    h = (double*) malloc(sizeof(double)*Nh);
    for(t=0;t<Nh;t++) scanf("%lf", &h[t]);

    printf("Type the input signal, size and data \n");
    scanf("%d", &Nx);
    x = (double*) malloc(sizeof(double)*Nx);
    for(t=0;t<Nx;t++) scanf("%lf", &x[t]);

    y = dspConv(x, Nx, h, Nh);

    printf("Result\n");
    vPrintf(stdout, y, Nx+Nh-1);
}




/* 
Calculates correlation between signals
uses convolution approach (algorithm not performance focussed)
*/

void
Correlation_Demo(){
    double *a, *b, *y; /* entrada, entrada, saida */
    int Na, Nb; /* tamanho dos sinais */
    int t; /* counters */

    printf("First signal: size and data \n");
    scanf("%d", &Na);
    a = (double*) malloc(sizeof(double)*Na);
    for(t=0;t<Na;t++) scanf("%lf", &a[t]);

    printf("Second signal: size and data\n");
    scanf("%d", &Nb);
    b = (double*) malloc(sizeof(double)*Nb);
    for(t=0;t<Nb;t++) scanf("%lf", &b[t]);

	/* uses convolution to calculate the correlation */
    y = dspCorr(a, Na, b, Nb);

    printf("Output signal \n");
    vPrintf(stdout, y, Na+Nb-1);

}

/*
Correlation Demos
same size signal correlation 
Uses mathlab code approach also not very effient (time performance)
*/

void
Correlation_MatlabCode_Demo(){
    double *a, *b, *y; /* entrada, entrada, saida */
    int N; /* tamanho dos sinais */
    int t; /* counters */

    printf("first signal: size and data\n");
    scanf("%d", &N);
    a = (double*) malloc(sizeof(double)*N);
    for(t=0;t<N;t++) scanf("%lf", &a[t]);

    printf("second signal: size and data\n");
    b = (double*) malloc(sizeof(double)*N);
    for(t=0;t<N;t++) scanf("%lf", &b[t]);

    /* -1 force the use of size maximum as number of lags */
    y = dspxCorr(a, b, N, -1);

    printf("output signal (xCorr Matlab)\n");
    vPrintf(stdout, y, 2*N-1);

    printf("output signal (dspCorr my C code)\n");
    y = dspCorr(a, N, b, N);
    vPrintf(stdout, y, 2*N-1);
}

/*
(You can use pipes < > to redirect the input/output from/to a file
Use the main to call any of the demo examples
*/

int main(){ 
   
	Dft_Fft_Demo();
	//Convolution_Demo();
	//Correlation_Demo();
	//Correlation_MatlabCode_Demo();
	printf("\n");
	system("PAUSE");

    return 0;
}


