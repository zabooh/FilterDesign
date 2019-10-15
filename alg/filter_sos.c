/* sos.c - IIR filtering by single second order section */

#include "filter_sos.h"

/* \(a, b, w\) are 3-dimensional */

/* \(a[0]=1\) always */


float FLT_float_sos(float *a, float *b, float *w, float x, int N) {
    float y;
    int i;

    for (i = 0; i < N; i++) {
        w[0] = x - a[1] * w[1] - a[2] * w[2];
        y = b[0] * w[0] + b[1] * w[1] + b[2] * w[2];

        w[2] = w[1];
        w[1] = w[0];

        a += 3;
        b += 3;
        w += 3;

        x = y;
    }

    return y;
}

long FLT_long_sos(long *a, long *b, long *w, long x, int N) {
    long y;
    int i;

    for (i = 0; i < N; i++) {
        w[0] = x - FRMULT32(a[1], w[1]) - FRMULT32(a[2], w[2]);
        y = FRMULT32(b[0], w[0]) + FRMULT32(b[1], w[1]) + FRMULT32(b[2], w[2]);

        w[2] = w[1];
        w[1] = w[0];

        a += 3;
        b += 3;
        w += 3;

        x = y;
    }

    return y;
}

