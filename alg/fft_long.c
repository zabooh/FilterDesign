/* fft.c  --  in-place decimation-in-time FFT */

//#include <plib.h>
#include <stdint.h>
#include <stdio.h>
#include "fft_long.h"


/* these arrays will be used by the FFT function */
long_complex FFT32BitInplace_Buffer[FFT_LENGTH]; /* complex input data */
uint32_t FFT32BitInplace_PowerSpectrum[FFT_LENGTH / 2]; /* power spectral density */
FFT_LONG_NUMERIC_TYPE FFT32BitInplace_StepResponse[FFT_LENGTH];

uint32_t FFT_long_GetLength(void) {
    return FFT_LENGTH;
}

void FFT_long_PutData(FFT_LONG_NUMERIC_TYPE Sample, uint32_t Index) {
    FFT32BitInplace_StepResponse[Index] = Sample; //(Sample >> (16 - Q31X_PRECISION));
    FFT32BitInplace_Buffer[Index].x = Sample; //(Sample >> (16 - Q31X_PRECISION));
    FFT32BitInplace_Buffer[Index].y = 0;
}

int32_t FTT_long_GetSampleData(uint32_t Index) {
    return FFT32BitInplace_StepResponse[Index];
}

void FFT_long_Calculate(void) {
    unsigned int Index;
    float ftemp;
    float t1, t2, t3, t4;

    //mPORTAToggleBits(BIT_14);
    FFT_long_fft(FFT_LENGTH, FFT32BitInplace_Buffer);
    //mPORTAToggleBits(BIT_14);

    /* compute the magnitudes */
    for (Index = 0; Index < FFT_LENGTH / 2; Index++) {
        /* compute the power spectral density of the FFT
         * using the sqrtf maths library function
         */
        t1 = FFT32BitInplace_Buffer[Index].x;
        t2 = FFT32BitInplace_Buffer[Index].y;
        t3 = t1 * t1;
        t4 = t2 * t2;
        ftemp = sqrtf(t3 + t4);
        if (ftemp == 0) {
            FFT32BitInplace_PowerSpectrum[Index] = 0;
        } else {
            FFT32BitInplace_PowerSpectrum[Index] = (uint32_t) (20.0 * log10f(ftemp));
        }
    }

}

uint32_t FFT_GetSpectrum(uint32_t Index) {
    return FFT32BitInplace_PowerSpectrum[Index];
}

void FFT_long_fft(int N, long_complex *X) {
    FFT_long_shuffle(N, X);
    FFT_long_dftmerge(N, X);
}

/* \(N\) must be a power of 2 */
void FFT_long_shuffle(int N, long_complex *X) {
    int n, r, B = 1;

    /* \(B\) = number of bits */
    while ((N >> B) > 0) {
        B++;
    }

    /* \(N = 2\sp{B}\) */
    B--;

    for (n = 0; n < N; n++) {
        /* bit-reversed version of \(n\) */
        r = FFT_long_bitrev(n, B);
        /* swap only half of the \(n\)s */
        if (r < n) continue;
        /* swap by addresses */
        FFT_long_swap(X + n, X + r);
    }
}

void FFT_long_dftmerge(int n, long_complex *XF) {
    double pi = 4. * atan(1.0);
    int k, i, p, q, m;
    long_complex A, B, V, W;

    m = 2;
    /* two \((M/2)\)-DFTs into one \(M\)-DFT  */
    while (m <= n) {
        /* order-\(M\) twiddle factor */
        //        mPORTAToggleBits(BIT_15);
        W = FFT_long_cmexp(FFT_long_cmplx(FFT_Q31X32(0.0), FFT_Q31X32(-2.0 * pi) / m));
        /* successive powers of \(W\) */
        V = FFT_long_cmplx(FFT_Q31X32(1.0), FFT_Q31X32(0.0));
        //        mPORTAToggleBits(BIT_15);
        /* index for an \((M/2)\)-DFT */
        for (k = 0; k < m / 2; k++) {
            /* \(i\)th butterfly; increment by \(M\) */
            for (i = 0; i < n; i += m) {
                /* absolute indices for */
                p = k + i;
                /* \(i\)th butterfly */
                q = p + m / 2;
                A = XF[p];
                /* \(V = W\sp{k}\) */
                B = FFT_long_cmul(XF[q], V);
                /* butterfly operations */
                XF[p] = FFT_long_cadd(A, B);
                XF[q] = FFT_long_csub(A, B);
            }
            /* \(V = VW = W\sp{k+1}\) */
            V = FFT_long_cmul(V, W);
        }
        /* next stage */
        m = 2 * m;
    }
}

int FFT_long_bitrev(int n, int B) {
    int m, r;

    for (r = 0, m = B - 1; m >= 0; m--)
        /* if \(2\sp{m}\) term is present, then */
        if ((n >> m) == 1) {
            /* add \(2\sp{B-1-m}\) to \(r\), and */
            r += TWO(B - 1 - m);
            /* subtract \(2\sp{m}\) from \(n\) */
            n -= TWO(m);
        }
    return (r);
}

void FFT_long_swap(long_complex *a, long_complex *b) {
    register long_complex t;

    t.x = (*a).x;
    t.y = (*a).y;
    (*a).x = (*b).x;
    (*a).y = (*b).y;
    (*b).x = t.x;
    (*b).y = t.y;

}

/* z = cmplx(x,y) = x+jy */
long_complex FFT_long_cmplx(FFT_LONG_NUMERIC_TYPE x, FFT_LONG_NUMERIC_TYPE y) {
    long_complex z;

    z.x = x;
    z.y = y;

    return z;
}

/* complex addition */
long_complex FFT_long_cadd(long_complex a, long_complex b) {
    return FFT_long_cmplx(a.x + b.x, a.y + b.y);
}

/* complex subtraction */
long_complex FFT_long_csub(long_complex a, long_complex b) {
    return FFT_long_cmplx(a.x - b.x, a.y - b.y);
}

/* complex multiplication */
long_complex FFT_long_cmul(long_complex a, long_complex b) {
    return
    FFT_long_cmplx(
            FFT_Q31X32_FRMULT(a.x, b.x) - FFT_Q31X32_FRMULT(a.y, b.y),
            FFT_Q31X32_FRMULT(a.x, b.y) + FFT_Q31X32_FRMULT(a.y, b.x)
            );
}

/* complex exponential */
long_complex FFT_long_cmexp(long_complex z) {
    float r = exp(FFT_Q31X32_TO_FLOAT(z.x));

    return
    FFT_long_cmplx(
            FFT_Q31X32(r * cos(FFT_Q31X32_TO_FLOAT(z.y))),
            FFT_Q31X32(r * sin(FFT_Q31X32_TO_FLOAT(z.y)))
            );
}

