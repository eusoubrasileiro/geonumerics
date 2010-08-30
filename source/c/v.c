/*
    V - LIBRARY
    Conjunto de funções de vetores especialmente
    desenvolvidas para utilização na solução de sistemas
    lineares pelo método de gaus com condensação pivotal
*/

#ifndef _V_H_
#define _V_H_


#include <stdio.h> /* fprintf ... */
/*
                            Funções para
                               vetores
*/

/*
imprime o vetor ...
*/
void
vPrintf(FILE *saida, double *v, int n){
    int i;
    fprintf(saida, "\n[ ");
    for(i=0; i<n; i++)
    fprintf(saida, "%3.2g ", v[i]);
    fprintf(saida, " ]\n");
}


/*
modulo simples de um numero
*/

double
vAbs(double valor){
    if(valor<0) valor*=-1;
    return valor;
}

/*
Multiplica todos os termos de um vetor (vetor) de tamanho (tamanho)
por um valor (valor)
*/

void
vMultv(double valor, int tamanho, double *vetor){
    register int i;
/* multiplica todos os termos da linha */
    for(i=0; i<tamanho; i++)
        vetor[i]*=valor;
}


/*
troca dois elementos de um vetor
dadas as posições dos elementos que deseja-se trocar
posorig e posdest
*/

void
vTrocav(double *vetor, int posorig, int posdest){
    double vaux; /* valor auxilira para a troca */

    /* copia o termo de origem em vaxuliar */
    vaux =  vetor[posorig];
    /* copia o de destino no de origem */
    vetor[posorig] = vetor[posdest];
    /* copia o auxialr em destino */
    vetor[posdest] = vetor[posorig];
}


/*
soma todos os termos de dois vetores/linhas
soma (vetora) com (vetorb) e coloca o resultado em (vsoma)
*/
void
vSomav(double *vetora, double *vetorb, double *vsoma, int vtamanho){
    register int j; /* percorre todos os termos dos vetores somando */
    /* soma todos os termos da linha */
    for(j=0; j<vtamanho; j++)
        vsoma[j] = vetora[j] + vetorb[j];
}

/*
copia um vetor em outro
*/

void
vCopiav(double *vetora, double *vetorcpy, int tamanho){
    register int j; /* percorre todos os termos do vetor copiando */
    /* copia ... */
    for(j=0; j<tamanho; j++)
        vetorcpy[j] = vetora[j];
}

/*
retorna indice do maior valor de um vetor
*/

int
vMaiorv(double *vetor, int tamanho){
    register int i; /* ora bolas um indice */
    int maior = 0; /* começa dizendo que o primeiro é o maior */

    for(i=1; i<tamanho; i++) /* se achar um elemento maior substitui */
        if(vetor[i]>vetor[maior]) maior = i;

    return maior;
}

/*
dado um vetor A e um vetor B
comprimentos nA e nB com nB = nA + deslocamento
desloca o vetor A de deslocamente e copia para B
*/

void
vDesloca(double *A, double*B, int nA, int deslocamento){
    int i;
    /*    deslocamento -> nA */
    for(i=0; i<nA; i++)
        B[i+deslocamento] = A[i];
}

/* ZERA */
void
vZera(double *v, int n){
    int i;
    for(i=0;i<n;i++) v[i]=0;
}

/*
make
f(t) = f(-t)
*/
#define SWAP(a,b) temp=(a);(a)=(b);(b)=temp

void
vTimeInvert(double *v, int n){
    int i, mid;
    double temp;

    /* central size of the vector*/
    mid=n/2;
    /* for vector size odd or even
     the number of swaps are n/2 */

    /* the last the first
      swap the first with the last
    the second with the last -1
     etc.. */
    if(n<1) /* no swapping */
        return;

    if(n==2){ /* just one swapping */
        SWAP(v[0],v[1]);
        return;
    }

    for(i=0; i<mid; i++)
        SWAP(v[i], v[n-1-i]);

}
#endif /* _V_H_ */
