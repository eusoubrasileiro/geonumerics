#include <stdio.h>
#include <stdlib.h>
#include "..\source\c\dsprocessing.c"


void
convolumain(int argc, char*argv[]){
    double *h, *x, *y; /* impulsiva, entrada, saida */
    int Nh, Nx; /* tamanho dos sinais */
    int t; /* counters */

    printf("resposta impulsiva do sistema, tamanho e dados\n");
    scanf("%d", &Nh);
    h = (double*) malloc(sizeof(double)*Nh);
    for(t=0;t<Nh;t++) scanf("%lf", &h[t]);

    printf("sinal de entrada do sistema, tamanho e dados\n");
    scanf("%d", &Nx);
    x = (double*) malloc(sizeof(double)*Nx);
    for(t=0;t<Nx;t++) scanf("%lf", &x[t]);

    y = pdsConv(x, Nx, h, Nh);

    printf("sinal de saida do sistema\n");
    vPrintf(stdout, y, Nx+Nh-1);
}

#define Nax 1000

void
dftmain(){
    unsigned int i=0, N;
    double A[Nax], an[Nax], bn[Nax];

    /* numero de amostras que vem, tudo o que estiver no arquivo */

    for(i=0; scanf("%lf", &A[i])==1; i++){} /* leh as amostras */
    N=i;
    fprintf(stderr,"n %u \n", N);
    
    /*    fprintf(stderr,"Namostras %lu\n", N); */
    pdsDft(N, A, an, bn); /* calcula os coeficientes */

    for(i=0; i<N; i++) /* imprime as belezinhas */
        printf("%lf %lf\n", an[i], bn[i]);

   //    for(i=0; i<N; i++) /* imprime as belezinhas */
   //        printf("%lf \n", sqrt(an[i]*an[i]+bn[i]*bn[i]));

}

/* any size signal */

void
corrmain(){
    double *a, *b, *y; /* entrada, entrada, saida */
    int Na, Nb; /* tamanho dos sinais */
    int t; /* counters */

    printf("primeiro sinal: tamanho e dados\n");
    scanf("%d", &Na);
    a = (double*) malloc(sizeof(double)*Na);
    for(t=0;t<Na;t++) scanf("%lf", &a[t]);

    printf("segundo sinal: tamanho e dados\n");
    scanf("%d", &Nb);
    b = (double*) malloc(sizeof(double)*Nb);
    for(t=0;t<Nb;t++) scanf("%lf", &b[t]);

    y = pdsCorr(a, Na, b, Nb);

    printf("sinal de saida \n");
    vPrintf(stdout, y, Na+Nb-1);

}

/* same size signal correlation */

void
xcorrAndCorrmain(){
    double *a, *b, *y; /* entrada, entrada, saida */
    int N; /* tamanho dos sinais */
    int t; /* counters */

    printf("primeiro sinal: tamanho e dados\n");
    scanf("%d", &N);
    a = (double*) malloc(sizeof(double)*N);
    for(t=0;t<N;t++) scanf("%lf", &a[t]);

    printf("segundo sinal: dados\n");
    b = (double*) malloc(sizeof(double)*N);
    for(t=0;t<N;t++) scanf("%lf", &b[t]);

    /* -1 force the use of size maximum as number of lags */
    y = pdsxCorr(a, b, N, -1);

    printf("sinal de saida xCorr\n");
    vPrintf(stdout, y, 2*N-1);

    printf("sinal de saida Correlate\n");
    y = pdsCorr(a, N, b, N);
    vPrintf(stdout, y, 2*N-1);
}

int main(int argc, char *argv[]){ /* eh soh escolher qual main vc quer rodar ... */
    //convolumain(argc, argv);
    xcorrAndCorrmain();
//    xcorrmain();
 //   dftmain();
system("PAUSE");

    return 0;
}


