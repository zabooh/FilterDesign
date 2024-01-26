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


/* fft.c  --  in-place decimation-in-time FFT */

#include <math.h>
#include "fft_float.h"
#include "complex_float.h"

void FFT_float_shuffle(int N, float_complex *X); /* \(N\) must be a power of 2 */
void FFT_float_dftmerge(int N, float_complex *XF);
void FFT_float_swap(float_complex *a, float_complex *b);
int FFT_float_bitrev(int n, int B);

void FFT_float_fft(int N, float_complex *X) {
    FFT_float_shuffle(N, X);
    FFT_float_dftmerge(N, X);
}

int FFT_float_bitrev(int n, int B) {
    int m, r;

    for (r = 0, m = B - 1; m >= 0; m--)
        if ((n >> m) == 1) { /* if \(2\sp{m}\) term is present, then */
            r += TWO(B - 1 - m); /* add \(2\sp{B-1-m}\) to \(r\), and */
            n -= TWO(m); /* subtract \(2\sp{m}\) from \(n\) */
        }

    return (r);
}

void FFT_float_swap(float_complex *a, float_complex *b) {
    register float_complex t;
    register float_complex c;
    register float_complex d;


    t.x = (*a).x;
    t.y = (*a).y;
    (*a).x = (*b).x;
    (*a).y = (*b).y;
    (*b).x = t.x;
    (*b).y = t.y;

    /*
       t = *a;
     *a = *b;
     *b =  t;
     */

}

void FFT_float_shuffle(int N, float_complex *X) /* \(N\) must be a power of 2 */ {
    int n, r, B = 1;

    while ((N >> B) > 0) /* \(B\) = number of bits */
        B++;

    B--; /* \(N = 2\sp{B}\) */

    for (n = 0; n < N; n++) {
        r = FFT_float_bitrev(n, B); /* bit-reversed version of \(n\) */
        if (r < n) continue; /* swap only half of the \(n\)s */
        FFT_float_swap(X + n, X + r); /* swap by addresses */
    }
}

void FFT_float_dftmerge(int N, float_complex *XF) {
    double pi = 4. * atan(1.0);
    int k, i, p, q, M;
    float_complex A, B, V, W;

    M = 2;
    while (M <= N) { /* two \((M/2)\)-DFTs into one \(M\)-DFT  */
        W = CPLX_float_cmexp(CPLX_float_cmplx(0.0, -2 * pi / M)); /* order-\(M\) twiddle factor */
        V = CPLX_float_cmplx(1., 0.); /* successive powers of \(W\) */
        for (k = 0; k < M / 2; k++) { /* index for an \((M/2)\)-DFT */
            for (i = 0; i < N; i += M) { /* \(i\)th butterfly; increment by \(M\) */
                p = k + i; /* absolute indices for */
                q = p + M / 2; /* \(i\)th butterfly */
                A = XF[p];
                B = CPLX_float_cmul(XF[q], V); /* \(V = W\sp{k}\) */
                XF[p] = CPLX_float_cadd(A, B); /* butterfly operations */
                XF[q] = CPLX_float_csub(A, B);
            }
            V = CPLX_float_cmul(V, W); /* \(V = VW = W\sp{k+1}\) */
        }
        M = 2 * M; /* next stage */
    }
}
