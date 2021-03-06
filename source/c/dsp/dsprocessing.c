/*

Biblioteca de processamento digital de sinais

*/

#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "dspVec.c" /* should change someway */
#include "dspFFT.c"
#include "dspFilters.c"

/*
a alegria da convolu��o
http://www.jhu.edu/~signals/convolve/index.html

convolu��o = calculo da resposta de um sistema linear
invariante no tempo � partir de um sinal de entrada

ex:
um sistema tem um resposta impulsiva : h[t] = [1, 2, 1, 3]
ou seja, para um impulso a resposta eh a impulsiva

exemplos para outros sinais
entrada -> saida
[1] -> [1, 2, 1, 3] resposta impulsiva
[0, 1] -> [0, 1, 2, 1, 3] invari�ncia no tempo
[2] -> 2*[1, 2, 1, 3] linearidade

[2,3] -> 2*[1, 2, 1, 3]+3*[0, 1, 2, 1, 3] = [2, 4, 2, 6]
                                          + [0, 3, 6, 3, 9]
                                          = [2, 7, 8, 9, 9]
o mesmo que (2+x)*(1+2x+x^2+3x^3) = 2+7x+8x^2+9x^3+9x^4

sinal de entrada no sistema  x[t]
resposta impulsiva do sistema h[t]
resposta do sistema para o sinal de entrada y[t] - resultado da convolucao

y[t] = soma k(0->Nx) x[k]*h[t-k]

*/


/*
*h, *x, *y  impulsiva, entrada, saida
o numero de amostras do sinal de saida eh Nh+(Nx-1)
FAZ
X[t](*)H[t] onde (*) eh o operador da convolucao

*/

double *
dspConv(double *x, unsigned int Nx, double *h, unsigned int Nh){
    unsigned int k; /* counters */
    double *h_k; /* ajudinha p/ a convolucao , guarda h[t-k]*/
    double *y; /* saida */
    unsigned int Ny;

    /* saida e help */
    Ny = Nh + (Nx-1); /* numero de amostras do sinal de saida */
    y = (double*) malloc(sizeof(double)*Ny);
    h_k = (double*) malloc(sizeof(double)*Ny);

    /* zera a saida first */
    vZera(y, Ny);

    /* a alegria da convolucao */
    for(k=0; k<Nx; k++){
        vZera(h_k, Ny); /* zera o vetor de ajuda */
/*        vPrintf(stdout, h_k, Ny); */
        vDesloca(h, h_k, Nh, k); /* copia o h[t] deslocado, ou seja, h[t-k] */
/*        vPrintf(stdout, h_k, Ny); */
        vMultv(x[k], Ny, h_k); /* muliplicada o h[t-k] por x[k] */
/*        vPrintf(stdout, h_k, Ny); */
        vSomav(y, h_k, y, Ny); /* soma com os valores anteriores */
/*        vPrintf(stdout, y, Ny); */
    }

    return y;
}

/*
Make the correlation between two traces or "vectors"


signal A a[t]
size : Na
signal B b[t]
size : Nb

correlation y[t]
size : Na+Nb-1

lag 0, in Nb

y[t] = soma k(-Nx->Nx) x[k]*h[t+k]

call convolution with h[-t], time inverted

Corr(a,b) != Corr(b,a)
because
Conv(a,b(-t)) != Conv(b, a(-t))

*/

double*
dspCorr(double *a, unsigned int sa, double *b, unsigned int sb){

    /* invert sempre a segunda */
    vTimeInvert(b, (int) sb);

    /* and call the convolution */
    /* a(t) * b(-t) */
    /* looks like that make a convolution of two vector of diferent size
    doesnt make sense because u would be dealing with a mess in the
    frequency domain, you first have to put than in the sample sample
    rate and after padded the smaller so they have the same size 
	Note: except for recursive filters? */
    return dspConv(a, sa, b, sb);

}

/*
Copy from mathlab
function xCorr, another way to do the same
thing above, but with arrays with the same size?? why the result is
diferent so?
TESTED with mathlab
N = size
m = delay
*/


double
_Rxy(double *x, double *y, int N, int m){
    int i;
    double res=0;
    /* m can displace max equal N-1 , N-m = 1*/
    for(i=0; i<N-m; i++)
        res += x[i+m]*y[i];

    return res;
}

double
_Ryx(double *y, double *x, int N, int m){
    int i;
    double res=0;
    /* m can displace max equal N-1 , N-m = 1*/
    for(i=0; i<N-m; i++)
        res += y[i+m]*x[i];

    return res;

}

double *
dspCorrx(double *x, double *y, int N, int m){
    int i;
	double *xcorr;

    /* if exceed bounds, put bounds
       if less, put bounds */
    m = ( m > N-1 || m < 0)? N-1 : m;

    xcorr = (double*) malloc(sizeof(double)*(2*m+1));

    for(i=-m; i<=m; i++) /* cover 2m+1 */
        xcorr[m+i] = (i>=0)? _Rxy(x, y, N, i): _Ryx(y, x, N, -i);

	return xcorr;
}

/*

dft - direct forier transform

*/



/*

Zn =

*/
/*
Calcula a transformada de fourier
de um vetor de dados
at� 4 bilh�es de dados 2^32 (unsigned int) = 4294967296
produto interno definido para duas fun��es complexas,
a e b, com N dados

<a(n), b(n)> = sumat�ria (n 0 -> N-1) a.b*
com b* conjugado de b

Considera-se a base de fun��es exponenciais, com N termos tal que:
exp[ (2pi/N)*k*t ]  k = 0, N-1 (s�o os vetores da base ortogonal)

2pi/N -> garante exponenciais com frequencias multiplas da frequencia
fundamental do sinal que eh assumida como 1/N. (Isso tambem � uma an�lise harm�nica
complexa)

Assim a proje��o � feita soh com a parte positiva
das exponencias pois a parte negativa � sim�trica.
Para o caso mais especifico do sinal de entrada real ainda:
o espectro � simetrico isso �
F(0->N/2) = F(N->N/2) ou = F(-N/2->N)
tambem
os coeficentes ncomplexs positivos s�o o
conjugado dos negativos, Z+n = (Z-n)*
o mesmo que F*(f) = F(-f)
|F(f)| funcao par
 Teta(f) fun��o �mpar

Tested with numpy FFT. Equal ABS values for rand noises. and PHASE
Totally validated!!

Ndiscrete fourier transform not normalized by N

also not very efficient due not taking the oportunity of having a simetrical spectrum

Should implement the inverse (ifft) but is meaningless
since I'm not going to use.

*/

#if BUILDING_DLL
# define EXPORTING __declspec (dllexport) __stdcall
#else /* Not BUILDING_DLL */
# define EXPORTING 
#endif /* Not BUILDING_DLL */

/*
good sense
*/

EXPORTING int
dspDft(double *samples, unsigned int ns, double *an, double *bn){
    unsigned int i, n; /* indice p/ somatoria e para os harmonicos */

    for(n=0; n < ns; n++){ /* 0 -> N-1 .. N frequencias */
        an[n]=bn[n]=0; /*to be used as accumulators * /
    /* produto interno de f(t) com cos(wn*i) e sin(wn*i) */
        for(i=0; i<ns; i++){ /* projeta nos cossenos e senos com t ou i indo de 0->N-1 */
            an[n] += samples[i]*cos( (double) n*2*M_PI*i/(ns));
            bn[n] += -1*samples[i]*sin( (double) n*2*M_PI*i/(ns));
        }
    }/* faz a projecao sobre o vetor exp[ (2pi/N)*k*t ] */
	return 0;
}

/*
Remember when calling to pass the adress of the pointer that will be alloced
double* an, bn;
dspDft_(samples, ns, &an, &bn);
*/
int
dspDft_(double *samples, unsigned int ns, double **an, double **bn){
   	double *as, *bs;    
    /* allocing to avoid problems */
	*an = (double*) malloc(sizeof(double)*ns);
	*bn = (double*) malloc(sizeof(double)*ns);
	as = an[0];
	bs = bn[0];   
	return dspDft(samples, ns, as, bs);
}

/*
imprime o espectro de frequencia em funcao da taxa de amostragem
e dos coefcientes complexox da transformada de fourier

tx de amostragem = espacamento em segundos entre as amostras

*/

void
dspEspectro(FILE *saida, double *a, double *b, unsigned int npontos, double txamos){
    unsigned int k;
    /* imprime */
    /* k -> cos(2pi*i*k/N) -> k = 0 frequencia k/(N*txamostragem) = 0 */
    /* k -> cos(2pi*i*k/N) -> k = 1 frequencia k/(N*txamostragem) (primeiro multiplo da frequencia fundamental) */

    /* coeficiente Epectro de frequencia */
    for(k=0; k<npontos; k++) /* precorre todos os coeficientes */
     fprintf(saida, "%.2f %.2f \n", k/(npontos*txamos), sqrt(a[k]*a[k]+b[k]*b[k]));
}

