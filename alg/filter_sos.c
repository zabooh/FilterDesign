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

