#ifndef __CMPLX_H
#define __CMPLX_H

#include "fft_long.h"

/* cmplx.h - complex arithmetic declarations */

//#include <math.h>                         /* in MSC and TC/BC, it declarares: */
/* \ttt{struct complex} and \ttt{cabs(z)} */
/* struct complex{float x, y;}; */ /* uncomment if neccessary */
/* float cabs(struct complex); */ /* uncomment if neccesary */

typedef struct {
    float x;
    float y;
} float_complex;


float_complex CPLX_float_cmplx(float x, float y); /* z = cmplx(x,y) = x+jy */
float_complex CPLX_float_conjg(float_complex z); /* complex conjugate of z=x+jy */
float_complex CPLX_float_cadd(float_complex a, float_complex b); /* complex addition */
float_complex CPLX_float_csub(float_complex a, float_complex b); /* complex subtraction */
float_complex CPLX_float_cmul(float_complex a, float_complex b); /* complex multiplication */
float_complex CPLX_float_rmul(float a, float_complex z); /* multiplication by real */
float_complex CPLX_float_cdiv(float_complex a, float_complex b); /* complex division */
float_complex CPLX_float_rdiv(float_complex z, float a); /* division by real */
float CPLX_float_real(float_complex z); /* real part Re(z) */
float CPLX_float_aimag(float_complex z); /* imaginary part Im(z) */
float_complex CPLX_float_cmexp(float_complex z); /* complex exponential */


#endif 
