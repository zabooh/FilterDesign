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

/* complex.c - complex arithmetic functions */

#include "complex_float.h"
#include <math.h>                         /* for MSC and TC/BC, it declares: */
/* \ttt{struct complex} and \ttt{cabs()} */
/* struct complex {float x, y;}; */ /* uncomment if not MSC or TC/BC */

/* uncomment if not MS or TC/BC */

/*  float cabs(z)
 *  complex z;
 *  {
 *      return sqrt(z.x * z.x + z.y * z.y);
 *  }
 */

float_complex CPLX_float_cmplx(float x, float y) /* z = cmplx(x,y) = x+jy */ {
    float_complex z;

    z.x = x;
    z.y = y;

    return z;
}

float_complex CPLX_float_conjg(float_complex z) /* complex conjugate of z=x+jy */
 {
    return CPLX_float_cmplx(z.x, -z.y); /* returns z* = x-jy */
}

float_complex CPLX_float_cadd(float_complex a, float_complex b) /* complex addition */ {
    return CPLX_float_cmplx(a.x + b.x, a.y + b.y);
}

float_complex CPLX_float_csub(float_complex a, float_complex b) /* complex subtraction */ {
    return CPLX_float_cmplx(a.x - b.x, a.y - b.y);
}

float_complex CPLX_float_cmul(float_complex a, float_complex b) /* complex multiplication */ {
    return CPLX_float_cmplx(a.x * b.x - a.y * b.y, a.x * b.y + a.y * b.x);
}

float_complex CPLX_float_rmul(float a, float_complex z) /* multiplication by real */ {
    return CPLX_float_cmplx(a * z.x, a * z.y);
}

float_complex CPLX_float_cdiv(float_complex a, float_complex b) /* complex division */ {
    float D = b.x * b.x + b.y * b.y;

    return CPLX_float_cmplx((a.x * b.x + a.y * b.y) / D, (a.y * b.x - a.x * b.y) / D);
}

float_complex CPLX_float_rdiv(float_complex z, float a) /* division by real */ {
    return CPLX_float_cmplx(z.x / a, z.y / a);
}

float CPLX_float_real(float_complex z) /* real part Re(z) */ {
    return z.x;
}

float CPLX_float_aimag(float_complex z) /* imaginary part Im(z) */ {
    return z.y;
}

float_complex CPLX_float_cmexp(float_complex z) /* complex exponential */ {
    float R = exp(z.x);

    return CPLX_float_cmplx(R * cos(z.y), R * sin(z.y));
}
