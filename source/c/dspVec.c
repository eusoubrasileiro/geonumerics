/*
    V - LIBRARY
    Conjunto de funções de vetores especialmente
    desenvolvidas para utilização na solução de sistemas
    lineares pelo método de gaus com condensação pivotal
*/

#ifndef _DSPVEC_H_
#define _DSPVEC_H_


#include <stdio.h> /* fprintf ... */
/*
                            Funções para
                               vetores
*/


/*

Simple Linear regression
For a set of pairs (x_i, y_i)
where x_i is given by the index
find the coeficients [a, b]
that the gives the minimum error
for
sum (yi - a - xi*b )**2 = E)
N number of pairs
return the coeficients a and b

b = sum xy - (sum  x * 1/N sum y )
b / = sum xx - [(sum x)^2 * 1/N ]
&
a = (1/N sum y)   - b * (1/N sum x)

*/

double* 
dspLinear(double* y, unsigned int N)
{
	double *a_b = (double*) malloc(sizeof(double)*2);
	unsigned int i;
	double xy, xx, x_, y_;
	a_b[0]=a_b[1]=xy=x_=y_=xx=0;

	for(i=0; i<N; i++){
		xy += i*y[i]; xx += i*i;
		x_ += i; y_ +=  y[i];
	}
	/*b = sum xy - (sum  x * 1/N sum y )
	b / = sum xx - [(sum x)^2 * 1/N ]*/
	a_b[1] = (xy - (x_*y_/N))/( xx - (x_*x_/N));
	a_b[0] = (y_/N) - (a_b[1]*x_/N);

	return a_b;
}


/*
remove the linear trend
from the data use the function above to
get a and b
make X = X - (a + b*xi)
where xi is the array index

also return the a_b coeficient pairs.
*/

double*
dspDetrendLinear(double* x, unsigned int N)
{
	unsigned int i;
	double* a_b = dspLinear(x, N);

	for(i =0; i<N; i++)
	{
		x[i] -= a_b[0] + a_b[1]*i;
	}
	
	return ab;
}

double*
dspRange(double begin, double end, double step)
{
	unsigned i, n = (unsigned int) ((end-begin)/step);
	
    if(end<begin || step==0)
        return NULL;
	
	// one sample more because its a closed interval
	double* range = new double[++n];

	for(i=0; i<n; i++)                
		range[i] = begin+i*step;

	return range;
}


/*
Sample the function passed, with values at array x
returns the array after sampling the function
*/
double*
dspSampleat(double *x, unsigned int n, double (*pfunc)(double x) )
{
    unsigned int i=0;
        /* closed interval one more sample */
    double *values = (double*) malloc(sizeof(double)*++n);
    for(i=0; i<n; i++)
        values[i] = pfunc(x[i]);
        
    return values;
}

/*
Sample the function passed.
Same above but with steps
from [x0, xf] with dx=step
*/

double*
dspSampleat_(double x0, double xf, double step, double (*func)(double x))
{
    unsigned int n = (unsigned int) (xf-x0)/step;
    double x;
    /* closed interval one more sample */
    double *values = (double*) malloc(sizeof(double)*++n);
    
    if(xf<x0 || step==0)
        return NULL;
    
    for(x=x0; x<=xf; x+=step)
        values[i] = pfunc(x);
    
    return values;        
}


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

	/* 
	example 5/2 = 2 
	swap(v[0],v[4])
	swap(v[1],v[3])
	example 6/2 = 3
	swap(v[0],v[5])
	swap(v[1],v[4])
	swap{v[2],v[3])
	*/
	for(i=0; i<mid; i++){ /* thousand of problems may arise it
						  it snot for this '{''}' guys */
        SWAP(v[i], v[n-1-i]);
	}

}
#endif /* _DSPVEC_H_ */
