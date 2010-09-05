
#ifndef _DSPROCESSING_H_
#define _DSPROCESSING_H_

extern double* dspLinear(double* y, unsigned int N);

extern void dspDetrendLinear(double* x, unsigned int N);

extern int dspDft(double *x, unsigned int ns, double *a, double *b);

extern int dspFft(double *x, unsigned int ns, double *a, double *b);

#endif /* _DSPROCESSING_H_ */
