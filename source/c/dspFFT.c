/*

Biblioteca de processamento digital de sinais
Fast Fourier Part

*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "dspVec.c" /* should change someway */

#ifndef _DSPFFT_H_
#define _DSPFFT_H_


#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif


typedef struct tagcomplex
{	
	double re, img;

} complex;


/*
imaginary exponential
return the complex equivalent
*/
complex
cExpImg(double c)
{
	complex ret;
	ret.re = cos(c);
	ret.img = sin(c);
	return ret;
}

complex
cMult(complex a, complex b)
{
	complex ret;
	ret.re = a.re*b.re-a.img*b.img;
	ret.img = a.re*b.img+a.img*b.re;
	return ret;
}

complex
cSum(complex a, complex b)
{
	complex ret;
	ret.re = a.re+b.re;
	ret.img = a.img+b.img;
	return ret;
}

/* could be a template for any array type */

complex*
cJoin(complex *Ca, unsigned int na, complex *Cb, unsigned int nb)
{
	unsigned int i;
	complex* ret = (complex*) malloc(sizeof(complex)*(na+nb));	

	for(i=0; i<na; i++)
		ret[i] = Ca[i];
	for(; i<na+nb; i++)
		ret[i] = Cb[i-na];

	return ret;
}

// just for ns%2==0 signals
double*
vGetEvenSamples(double* Sx, unsigned int n)
{
	unsigned int i;
	double *even = (double*) malloc(sizeof(double)*n/2);
	for(i=0; i<n/2; i++) /* gets ns=5: indexes 0, 2, 4, 6, 8 .. */
		even[i] = Sx[i*2]; 
	return even;
}

// just for ns%2==0 signals
double*
vGetOddSamples(double* Sx, unsigned int n)
{
	unsigned int i;
	double *odd = (double*) malloc(sizeof(double)*n/2);
	for(i=0; i<n/2; i++) /* gets ns=5: 1, 3, 5, 7, 9 .. */
		odd[i] = Sx[i*2+1]; 
	return odd;
}

/*
-1 or 1 for true or false
*/

int 
isPower2(unsigned int x)
{
    /*
     x And x-1       
     8 = 01000   (x)
     7 = 00111  (x-1)
     & = 00000   logo é power 2     
    */
	return ( (x > 0) && ((x & (x - 1)) == 0) ) ? 1 : -1;
}

/*
what u expect by the name
*/

unsigned int
nextPower2(unsigned int n)
{
	unsigned int powof2 = 1; /* start with 2**0 = 1, and shift bits to the left */ 

	while( powof2 < n ) 
		powof2 <<= 1;

	return powof2;
}



/*
Fft coley Tukey recursive way
just for real data input,
optimzed to simetrical spectrum ( automatic from Coley tukey algorithm definiton)
Just power two size's allowed to enter here
*/

complex*
dspInternFft(double* samples, unsigned int ns)
{	
	unsigned int k;
	complex *X, *Even, *Odd;
	double *odd, *even;
	complex _minus;
	/* just to be used to multiply for -1 */
	_minus.re = -1;
	_minus.img = 0;
	/* the fft result complex and imaginary parts */
	
	if(ns == 1) /* simplest case */
	{
		// just here we need to alloca in the other case we already have
		// it allocated by the Join Funct
		X = (complex*) malloc(sizeof(complex));
		X[0].re = samples[0]; /* samples[0]*exp(??? * 0) = samples[0] */
		X[0].img = 0;
	}
	else
	if(ns == 2) /* 2nd simplest case */
	{
		X = (complex*) malloc(sizeof(complex)*ns);
		X[0].re = samples[0] + samples[1];
		X[1].re = samples[0] - samples[1]; 
		X[0].img = X[1].img = 0;

		return X;
	}  
    else
	{
		even = vGetEvenSamples(samples, ns);
		odd = vGetOddSamples(samples, ns);
		Even = (complex*) dspInternFft(even, ns/2);
		Odd = (complex*) dspInternFft(odd, ns/2);
		X = (complex*) malloc(sizeof(complex)*ns);

		for(k = 0 ; k<ns/2; k++){
            /* Even : even part */ 
			/* Odd : odd part */
			/* could improve performance here removing all those function calls to simple math */
            X[k] = cSum(Even[k], cMult( cExpImg((double) -1*2*M_PI*k/ns), Odd[k]) );
            X[k+(ns/2)] = cSum(Even[k], cMult( cMult( cExpImg((double) -1*2*M_PI*k/ns), Odd[k]) , _minus));
		}
		free(even);
        free(odd);		
		free(Even);
		free(Odd);
	}
	/*
    This bellow is already exploited above so no more simetry to explore
    create the other half N/2 -> N, just simetric copy around ns_orig/2 
	
	for input signal size N = 4	
	X[3]=*X[2-1]; (index 2 point of simetry)
	
	where * is conjugate o X

	for input signal size N=8
	X[5]=*X[4-1]; (index 4 point of simetry)
	X[6]=*X[4-2]; (index 4 point of simetry)
	X[7]=*X[4-3]; (index 4 point of simetry)
    */

	return X;
}


/*
Callable function for real input signals
must be power of 2 length
*/

/*
good sense
*/

int
dspFft(double *samples, unsigned int ns, double *an, double *bn)
{
	unsigned int i; 
	complex* ret;

	if(isPower2(ns)==-1)
		return -1;

	/* already power 2 just make the math */
	ret = dspInternFft(samples, ns);

	/* copy back in the coeficient form just the original size */
	for(i=0; i<ns; i++){
		an[i] = ret[i].re;
		bn[i] = ret[i].img;
	}

	free(ret);
	return 0;
}

/*
Remember when calling to pass the adress of the pointer that will be alloced
double* an, bn;
dspFft_(samples, ns, &an, &bn);
*/
int
dspFft_(double *samples, unsigned int ns, double **an, double **bn){    
   	double *as, *bs;
    /* allocing to avoid problems */
	*an = (double*) malloc(sizeof(double)*ns);
	*bn = (double*) malloc(sizeof(double)*ns);
	as = an[0];
	bs = bn[0]; 
	return dspFft(samples, ns, as, bs);
}

#endif /* _DSPFFT_H_ */
