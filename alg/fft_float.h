#ifndef __FFT_H
#define __FFT_H

/* fft.c  --  in-place decimation-in-time FFT */

#include "complex_float.h"

void FFT_float_fft(int N, float_complex *X);

#endif
