/*

Biblioteca de processamento digital de sinais

*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "dsprocessing.h"
#include "v.c" /* should change someway */

/*
a alegria da convolução
http://www.jhu.edu/~signals/convolve/index.html

convolução = calculo da resposta de um sistema linear
invariante no tempo à partir de um sinal de entrada

ex:
um sistema tem um resposta impulsiva : h[t] = [1, 2, 1, 3]
ou seja, para um impulso a resposta eh a impulsiva

exemplos para outros sinais
entrada -> saida
[1] -> [1, 2, 1, 3] resposta impulsiva
[0, 1] -> [0, 1, 2, 1, 3] invariância no tempo
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

double *
dspCorr(double *a, unsigned int sa, double *b, unsigned int sb){

    /* vPrintf(stdout, b, sb); */ /* pra checar se o time invert ta certo e ta */
    /* invert sempre a segunda */
	/* not sure if this time invert is right */
    vTimeInvert(b, (int) sb);
    //vPrintf(stdout, b, sb);

    /* and call the convolution */
    /* a(t) * b(-t) */
    /* looks like that make a convolution of two vector of diferent size
    doesnt make sense because u would be dealing with a mess in the
    frequency domain, you first have to put than in the sample sample
    rate and after padded the smaller so they have the same size */
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
Rxy(double *x, double *y, int N, int m){
    int i;
    double res=0;
    /* m can displace max equal N-1 , N-m = 1*/
    for(i=0; i<N-m; i++)
        res += x[i+m]*y[i];

    return res;
}

double
Ryx(double *y, double *x, int N, int m){
    int i;
    double res=0;
    /* m can displace max equal N-1 , N-m = 1*/
    for(i=0; i<N-m; i++)
        res += y[i+m]*x[i];

    return res;

}

/*
[2, 3] e [1, 2, 3, 4, 1]
na correlacao normal seria




*/
double *
dspxCorr(double *x, double *y, int N, int m){
    int i;
	double *xcorr;

    /* if exceed bounds, put bounds
       if less, put bounds */
    m = ( m > N-1 || m < 0)? N-1 : m;

    xcorr = (double*) malloc(sizeof(double)*(2*m+1));

    for(i=-m; i<=m; i++) /* cover 2m+1 */
        xcorr[m+i] = (i>=0)? Rxy(x, y, N, i): Ryx(y, x, N, -i);

	return xcorr;
}

/*

dft - direct forier transform

*/

#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif

/*

Zn =

*/
/*
Calcula a transformada de fourier
de um vetor de dados
até 4 bilhões de dados 2^32 (unsigned int) = 4294967296
produto interno definido para duas funções complexas,
a e b, com N dados

<a(n), b(n)> = sumatória (n 0 -> N-1) a.b*
com b* conjugado de b

Considera-se a base de funções exponenciais, com N termos tal que:
exp[ (2pi/N)*k*t ]  k = 0, N-1 (são os vetores da base ortogonal)

2pi/N -> garante exponenciais com frequencias multiplas da frequencia
fundamental do sinal que eh assumida como 1/N. (Isso é uma análise harmônica
complexa)

Assim a projeção é feita soh com a parte positiva
das exponencias pois espera-se, que como para toda função real
o espectro é par e a fase ímpar não precisamos fazer a projeção
de -inf -> +inf ou de -N+1 -> N-1 basta a parte positiva (ou negativa)

os coeficentes ncomplexs positivos são o
conjugado dos negativos, Z+n = (Z-n)*
o mesmo que F*(f) = F(-f)
|F(f)| funcao par
 Teta(f) função ímpar

Tested with numpy FFT. Equal ABS values for rand noises. and PHASE
Totally validated!!
*/

void
dspDft(unsigned int Namostras, double *Amostras, double *an, double *bn){
    unsigned int i, n; /* indice p/ somatoria e para os harmonicos */


    for(n=0; n < Namostras; n++){ /* 0 -> N-1 .. N frequencias */
        an[n]=bn[n]=0; /*to be used as accumulators * /
    /* produto interno de f(t) com cos(wn*i) e sin(wn*i) */
        for(i=0; i<Namostras; i++){ /* projeta nos cossenos e senos com t ou i indo de 0->N-1 */
            an[n] += Amostras[i]*cos( (double) n*2*M_PI*i/(Namostras));
            bn[n] += -Amostras[i]*sin( (double) n*2*M_PI*i/(Namostras));
        }
    }/* faz a projecao sobre o vetor exp[ (2pi/N)*k*t ] */
}

/*
inversa..

*/

void
dspIDft(unsigned int Namostras, double *Amostras, double *an, double *bn){
    unsigned int i, n; /* indice p/ somatoria e para os harmonicos */

    for(n=0; n < Namostras; n++){ /* 0 -> N-1 .. N frequencias */
        an[n]=bn[n]=0; /*to be used as accumulators * /
    /* produto interno de f(t) com cos(wn*i) e sin(wn*i) */
        for(i=0; i< Namostras; i++){ /* projeta nos cossenos e senos com t ou i indo de 0->N-1 */
            an[n] += Amostras[i]*cos( (double) n*2*M_PI*i/(Namostras));
            bn[n] += Amostras[i]*sin( (double) n*2*M_PI*i/(Namostras));
        }
    }/* faz a projecao sobre o vetor exp[ (2pi/N)*k*t ] */
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

double* dspLinear(double* y, unsigned int N)
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
*/

void dspDetrendLinear(double* x, unsigned int N)
{
	unsigned int i;
	double* a_b = dspLinear(x, N);

	for(i =0; i<N; i++)
	{
		x[i] -= a_b[0] + a_b[1]*i;
	}
}



/*
IIR - infinite impulse response

Monta um filtro "caixa" no tempo, dados
1) A frequencia de corte Fc
2) O numero de amostras do filtro
3) A taxa de amostragem do filtro (dt) (deve ser igual ao do sinal à ser filtrado)
pq senao vc tah mulptiplicando outra coisa na frequencia...

A expressao para o filtro é (inversa da caixa na frequencia)

Box(t)  =  ( A sin(2 pi Fc t) ) / ( pi t )

A = dt (tx de amostragem)

O tamanho do filtro no tempo eh ( N dt )
A amostragem do filtro no tempo é feito de t:
-(N dt)/2 ate (N dt)/2
para garantir a simetria do filtro nao modificando seu espectro
o Box(0) no limite eh igual à

*/

double
*dspIIRBox(double Fc, unsigned int Na, double dt){

    return 0;
}

