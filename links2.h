#ifndef __LINKSJG_H_
#define __LINKSJG_H_

#include <stdbool.h>

double PI();
void  fb(float a, float b, float c, float d, float T1, float T2, float *T3, float *T4, bool crossed);
void cs(float a, float b, float c, float T1, float T2, float *T3, float *d );
void ics(float a, float c, float d, float T2, float G, float *T3, float *T4, float *b );
void pAS(float O2x, float O2y, float T2, float d2, float s, float *Sx, float *Sy);
void pBP(float Ax, float Ay, float T3, float d3, float p, float *Px, float *Py);
void pBU(float O4x, float O4y, float T4, float d4, float u, float *Ux, float *Uy);
void tran(float T3, float T4, float *mu);
int grash(float a, float b, float c, float d);
void extran(float a, float b, float c, float d, float *mu1, float *mu2);
void tog(float a, float b, float c, float d, float *Ttog1, float *Ttog2);
void fbv(float a, float b, float c, float T2, float T3, float T4, float w2, float *w3, float *w4 );
void csv(float a, float b, float c, float T1, float T2, float T3, float w2, float *w3, float *ddt);
void csa(float a, float b, float c, float T1, float T2, float T3, float w2, float w3, float al2, float *al3, float *ddtt );
void fba(float a, float b, float c, float T2, float T3, float T4, float w2, float w3, float w4, float al2, float *al3, float *al4 );





#endif
