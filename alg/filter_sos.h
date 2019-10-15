#ifndef _SOS_H
#define _SOS_H

#include "..\FilterDesign\filter_design.h"

/* sos.c - IIR filtering by single second order section */

#define FRMULT32(a,b)    ((long)(((long long)a * (long long)b)>>INTEGER_PRECISION))
#define QX(X)   ((X < 0.0) ? (long)((1<<INTEGER_PRECISION)*(X) - 0.5) : (long)((1<<INTEGER_PRECISION)*(X) + 0.5)) 
#define QX_ONE   (1<<INTEGER_PRECISION)


float FLT_float_sos(float *a, float *b, float *w, float x, int N);
long FLT_long_sos(long *a, long *b, long *w, long x, int N);

#endif
