#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "dspVec.c" /* should change someway */
#include "dspFFT.c" /* FFT and IFFT */

#ifndef _DSPVEC_H_
#define _DSPVEC_H_


#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif

/* filter convolution options */
#define DETREND_LINEAR 1
#define HANNING_WINDOW 1

// <summary>
//"""
//Gives the number of samples (odd number) needed to sample our filter operator
//based on:
//1) the relative transition bandwidth RTbtw WE WANT 
//(Transition band size divide by cut-off frequency)
//2) Filter cut-off frequency 
//3) Sample rate that's gonna be used to sample the filter
//This aproximation is valid for Sinc Filters.
//a RTbtw = 20% is a reasonable value for a non distorded frequency response	
//"""

/* 
gets the ideal filter size for the desired relative bandwidth transition 
unsigned int (size)
rtbtw desired relative bandwidth transition (0.2 is normally good)
fc cuttof frequency (filter frequency cuttof)
dt sample rate for the filter
*/
#define filter_size(rtbtw, fc, dt) ((unsigned int) (2/(rtbtw)*(fc)*(dt))
/*
gets the next odd integer 
unsigned int input
*/
#define next_odd(n) ((n)%2==0)? (n)+1 : (n);




double*
dspFilterConvFft(double* signal, unsigned int ns, 
    double* filter_kernel, unsigned int nf, 
    int detrend_option)
{
	unsigned int i, nf_2; /* counter */
    double* linear_ab; /* original trend removed from the signal (coeficients a+b*x */	
    
    if(detrend_option==DETREND_LINEAR){    
    	/* first remove the stationary to
        avoid problems that might exists, linear removal */
    	linear_ab = dspDetrendLinear(signal, ns);
    	/* store the original trend */    
    }
    
    /* filter kernel must be odd to have a simetric convolution */
    if(nf%2 == 0)    
        return null;    
        
    Sor = ns;

    if(ns%2 == 0)//must be odd
        return NULL;

    nf_2 = (nf-1)/2; // half size of the filter

    // easy part just padding zeros, with the half size that is bigger
    signal = dspAppend(signal, ns, dspZeros(nf_2+1), nf_2+1); 
    
    Ssize = ns+nf_2+1; // signal size, now signal is bigger, more samples

    // wrap around form modification
    // the above creates simetrecally around 0
    fir_hannEnd = fir_hann[0:Fhsize]; #size ...
    fir_hannBeg = fir_hann[Fhsize:]; # size ...
    // make the same size the input signal, considering  the signal size greater than the filter
    // case the signal is smaller than filter
    // putting zeros in the hole between the initial and end part of the wrap around form!
    if(Ssize-(Fhsize*2+1) >= 0):
        fir_hannBeg = numpy.append(fir_hannBeg, numpy.zeros(Ssize-(Fhsize*2+1)));

    #put in the wrap around form
    fir_hann = numpy.append(fir_hannBeg, fir_hannEnd);
    Sfsize = numpy.size(fir_hann);

    if(Ssize < Sfsize): # the else condition of above
        # they dont have the same suze
        # padd the input signal with zeros until they have the same size
        # to be able to perform the convolution
        signal = numpy.append(signal, numpy.zeros(Sfsize-Ssize));

    print "F End %d" % numpy.size(fir_hannEnd)
    print "F Beg %d" % numpy.size(fir_hannBeg)
    print "F %d" % numpy.size(fir_hann)
    print "S %d" % numpy.size(signal)

    sigrs = pylab.real(pylab.ifft(pylab.fft(signal)*pylab.fft(fir_hann)));

    #python is nice hehhe
    # signal = signal[:-(Fhsize+1)] #remove the last part added by the filter avoiding
    sigrs_ = sigrs[:Sor]; # get the just the first part equivalent to the signal size
    sigrs_ += trend; #add the trend back
    
    if(detrend_option==DETREND_LINEAR){    
    	/* 
        add the trend back to the signal (putting it as the orignal one)
        a + b*i
        */
    	for(i=0; i<ns; i++)
    	   signal += linear_ab[0] + linear_ab[1]*i;
    }
    
    return sigrs;
    
}

/*
    Hanning window working for smooth the frequency response
    removing ripples most used one
*/
double*
dspWindowHann(unsigned int n)
{    
    /* thats because that func returns a close interval */    
    return dspSampleat_(0, n-1, 1, double (*Fx)(double x){ return sin(x*M_PI/(n-1)) });
}

/*
filtro trapezoidal passa baixa ( rampa = Ramp, Frequencia de corte = Fc )
infinite impulse response cortada ficando finita (fir) e filtro nao causal
1) filtro eh truncado no tempo pois nao podemos ter infinitas amostras...
logo o filtro eh multiplicado por uma funcao caixa no tempo que eh o
mesmo que convolver com o sinc dessa caixa na frequencia
que acarreta gibs no modulo e fase
2) o filtro eh amostrado no tempo, ou seja, multiplicado por deltas de
amplitude 1 e espacados de dt. Que eh o mesmo que convolver na frequencia
com deltas espacados de 1/dt e amplitude 1/dt
dt = dt  Deve ser  tal que  === 1/2*dt > Fc
*/

double*
dspIIR_SincTrapezLowPass(unsigned int n, double dt, double fc, double ramp, int window_option)
{
    unsigned int i;
	// this filter is the result of convolution of two box in frequency
	// the box A with its side/2 = a
	double a = Ramp/2;
	// th box B with its side/2 = b
	double b = fc - (ramp/2);

	// this doesnt work, ramp = 0, and other errors 
	/* "Impossible size of ramp" */
	if(a==0 || ramp > fc || fc < 1/dt)
	   return NULL;		

	double* x = NULL; // imput x/t values to sample the filter operator 
	double* filter = (double*) malloc(sizeof(double)*n); // filter kernel sampled
	double* window_taper;

	// Amostra simetricamento o operador do filtro em torno do zero
	// nao importa se o numero de pontos for par ou impar, exitem opcoes pros dois casos	
	
	// caso seja impar amostra incluindo o 0,  amostra perfeitamente simetrico em torno do zero
	// um workaround eh utilizado para evitar divsao por 0 e utiliza-se o limite em 0 para setar o valor no 0
	if (N % 2 != 0)
	{		
        /* simetric around 0 */
		x = dspRange(-dt*(n-1)/2, dt*(n-1)/2, dt);
		// para evitar excessao de divisao por 0		
		x[n / 2] = 1.0;
		// sample the function
		for (i = 0; i < n; i++)
			filter[i] = dt * sin(2 * M_PI * a * x[i]) * sin(2 * M_PI * b * x[i]) / (M_PI * M_PI * x[i] * x[i] * 2 * a);

		// inverso da convolucao de
		// duas caixas na frequencia
		// dt multiplicando serve para garantir a o espc. amplitude em 1	
		// a divisao por 2*a, serve para garantir o spec. amplitude em 1, 
		// pois o resultado da convolucao de duas caixas de amplitude 1 na frequencia [-a, a], [-b, b] com b > a 
		// resulta no trapezio com amplitude maxima igual 2*a		
		// set o valor de y no 0, baseado no limite x->0 do operador
		// fazendo o limite (derivando em baixo e em cima) chega-se para a => cos(2*pi*a*t)*2*a
		// com t =0 => 2*a, ou seja para a e b => 2*a*2*b
		// limite da multplicacao eh a multiplicacao dos limites
		filter[n / 2] = 2 * a * 2 * b * dt / (2 * a);
	}
	else
	{
		//amostra  simetricamente em torno do zero sem passar pelo 0, mais facil, nao precisa evitar divisao por zero		        
		x = dspRange(-dt * (double)(n - 1) / 2, dt * (double)(n - 1) / 2, dt);
		// sample the function
		for (i = 0; i < n; i++)
			filter[i] = dt * sin(2 * M_PI * a * x[i]) * sin(2 * M_PI * b * x[i]) / (M_PI * M_PI * x[i] * x[i] * 2 * a);
	}
	
	/* apply the taper in the filter kernel */
    if(window_option==HANNING_WINDOW)
    {
        /* gets the taper */
        window_taper = dspWindowHann(n);        
        for(i=0; i < n; i++) /* applies the taper */
            filter[i] *= window_taper[i];
        
    }

	// uhuu the filter done
	return filter;
}





double* 
dspIIR_SincTrapezLowPassApply(double *signal, unsigned int ns, double dt, double fc, double ramp)
{   
	// filter operator in time
	double* filter_kernel;
	unsigned int nf; /* number of samples for the filter*/
	
	/* Filter definition */
	// Size must be always ODD!! ALWAYS!! 2n+1 where n is the order of the filter            
	// Previous problem... if the sample rate is for example 1E-5 using 
	// a filter with 101 samples, cant represent the spectrum properly            
	// and transition bandwith, the length of the transition part,
	// for a sinc filter the transition btw can be aproximated by
	// BTW = (4/Nf)*(1/(2*dt)) 
	// e.g dt=2E-5 Nf=101 
	// BTW = 715 Hz  considering
	// that I want to filter 125 hz... ahhaa Horrible!! u dont even have a filter!
	// relative band width transition achived by dividing by the cutoff frequency            
	// 20% has an accepted ripple/transition
	// using also the hanning window!!
	nf = next_odd(filter_size(0.2, fc, dt));
	nf = (nf < 101) ? 101 : nf; // minimum size 101 ??      
	
	// Creates the filter Kernel
	filter_kernel = TrapezoidalLowPass(nf, dt, fc, ramp, HANNING_WINDOW);
	
	// applys the filter fft convolution clipping and whatever... doesnt alter
	// the signal input array
	signal = dspFilterConvFft(signal, ns, filter_kernel, nf, DETREND_LINEAR);

    return signal;	
}

/*
IIR - infinite impulse response

Monta um filtro "caixa" no tempo, dados
1) A frequencia de corte Fc
2) Ordem do filtro: On*2+1 = N:numero de amostras do filtro
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
*dspIIRBox(double Fc, unsigned int On, double dt){

    return 0;
}

#endif /* _DSPFILTERS_H_ */
