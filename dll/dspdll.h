#ifndef _DLL_H_
#define _DLL_H_

#if BUILDING_DLL
# define DLLIMPORT __declspec (dllexport) __stdcall
#else /* Not BUILDING_DLL */
# define DLLIMPORT __declspec (dllimport) __stdcall
#endif /* Not BUILDING_DLL */

#ifdef _MSC_VER
extern "C" {
#endif

DLLIMPORT 
double * pdsCorr(double *a, unsigned int sa, double *b, unsigned int sb);

DLLIMPORT 
double * pdsxCorr(double *x, double *y, int N, int m);

DLLIMPORT 
void pdsDft(unsigned int Namostras, double *Amostras, double *an, double *bn);


#ifdef _MSC_VER
}
#endif


#endif /* _DLL_H_ */
