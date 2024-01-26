/*********************************************************************

Copyright (c) 2019 Martin Ruppert

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

*********************************************************************/

#ifndef __FFT32_H
#define __FFT32_H

#include <stdint.h>
#include <math.h>

//#define __USE_FLOAT

#define FFT_LENGTH 512

#ifdef __USE_FLOAT
#define FFT_LONG_NUMERIC_TYPE	float
#define FFT_Q31X_FRMULT(a,b) a * b
#define FFT_Q31X32(x) x
#define FFT_Q31X32_TO_FLOAT(y) y
#define FFT_Q31X32_PRECISION 0
#else
#define FFT_LONG_NUMERIC_TYPE	int32_t
#define FFT_Q31X32_PRECISION 2
#define FFT_Q31X32_DECIMAL_POINT_POSITION ((32-2) - FFT_Q31X32_PRECISION)
#define FFT_Q31X32_FRMULT(a,b) \
    ( \
        ( \
            (int32_t)( \
            ((int64_t)a * (int64_t)b)>>32) \
        ) << (2 + FFT_Q31X32_PRECISION) \
    )
#define FFT_Q31X32(X) \
    ( \
        (X < 0.0) ? \
            (int32_t) ((1 << FFT_Q31X32_DECIMAL_POINT_POSITION)*(X) - 0.5) \
        : \
            (int32_t) ((1 << FFT_Q31X32_DECIMAL_POINT_POSITION)*(X) + 0.5) \
    )

#define FFT_Q31X32_TO_FLOAT(y)  ((float)y / (float)(1 << FFT_Q31X32_DECIMAL_POINT_POSITION))
#endif


#define TWO(x)  (1 << (x))

typedef struct {
    FFT_LONG_NUMERIC_TYPE x;
    FFT_LONG_NUMERIC_TYPE y;
} long_complex;

uint32_t FFT_long_GetLength(void);
void FFT_long_PutData(FFT_LONG_NUMERIC_TYPE Sample, uint32_t Index);
void FFT_long_Calculate(void);
int32_t FTT_long_GetSampleData(uint32_t Index);

void FFT_long_fft(int N, long_complex *X);
void FFT_long_shuffle(int N, long_complex *X); /* \(N\) must be a power of 2 */
void FFT_long_dftmerge(int N, long_complex *XF);
int FFT_long_bitrev(int n, int B);

void FFT_long_swap(long_complex *a, long_complex *b);
long_complex FFT_long_cmplx(FFT_LONG_NUMERIC_TYPE x, FFT_LONG_NUMERIC_TYPE y); /* z = cmplx(x,y) = x+jy */
long_complex FFT_long_cadd(long_complex a, long_complex b); /* complex addition */
long_complex FFT_long_csub(long_complex a, long_complex b); /* complex subtraction */
long_complex FFT_long_cmul(long_complex a, long_complex b); /* complex multiplication */
long_complex FFT_long_cmexp(long_complex z); /* complex exponential */


#endif
