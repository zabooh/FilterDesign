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

#ifndef __FILTER_DESIGN
#define __FILTER_DESIGN

#include "math.h"

#define bool unsigned int
#define false 0
#define true  1

#define arsinh(x)           (log((x)+sqrt((x)*(x)+1)))
#define arcosh(x)           (log((x)+sqrt((x)*(x)-1)))	// ACHTUNG nur für x-Werte >= 1 !!!
#define MAX                 16
#define MAXDIGFILTKOEFF     5*MAX
#define	DB                  log(10)/10
#define	PI                  3.14159265358979
#define	FALSE               0
#define	TRUE                1

#define INTEGER_PRECISION   27

/* Codes for the selection of the filter characteristic */
#define BUTT        0x01
#define TSCHE       0x02
#define CAU         0x04

/* Codes for selecting the filter type */
#define TP          0x08
#define HP          0x10
#define BP          0x20
#define BS          0x40

#define N_EQUAL_ZERO	1

typedef struct {
    int Error;
    float Samplerate;
    short Coefficient_Block_2[64];
    short Coefficient_Block_1[16];
    float Az_Biquad[MAX + 1];
    float Ay_Biquad[MAX + 1];
    float Ax_Biquad[MAX + 1];
    float Aw_Biquad[MAX + 1];
    float Av_Biquad[MAX + 1];
    float F_Coeff[MAX + 1];
    float E_Coeff[MAX + 1];
    float D_Coeff[MAX + 1];
    float C_Coeff[MAX + 1];
    float B_Coeff[MAX + 1];
    float A_Coeff[MAX + 1];
    float F4_Analog;
    float F3_Analog;
    float F2_Analog;
    float F1_Analog;
    float f4_Digital;
    float f3_Digital;
    float f2_Digital;
    float f1_Digital;
    float as_Blocking_Attenuation;
    float ad_Passband_Attenuation;
    float Epsilon;
    int Sub_Filter;
    int n_Order;
    float constant_values[20];
    float omega_d[20];
    float s_imaginary_part_complex_zeros[20];
    float r_real_part_complex_zeros[20];
    float omega_s;
    float cauer_const;
    float b_reference_lp_coeff[MAX + 1];
    float a_reference_lp_coeff[MAX + 1];
} DSP_DATA;



void FLD_Frequency_Transformation();
void FLD_Digitale_Koeffizienten_ermitteln(void);
void FLD_Coefficient_Assignment(void);
int FLD_BS_Transformation(int filt_char);
int FLD_BP_Transformation(int filt_char);
int FLD_Referenz_TP_Cauer();
int FLD_BS_Cauer();
int FLD_BP_Cauer();
int FLD_HP_Cauer();
int FLD_TP_Cauer();
void FLD_Get_Constants();
float FLD_gammaM(float a, float b);
float FLD_u_frq(float omega_s, int nn);
float FLD_u_db(float as, float ad);
float FLD_Sn(float u, float z);
void FLD_InitCoeffsFloat(float *num, float *den, float *taps, int N);
void FLD_InitCoeffsFixpoint(long *num, long *den, long *taps, int N);
int FLD_Referenz_TP_Tschebycheff();
int FLD_BS_Tschebycheff();
int FLD_BP_Tschebycheff();
int FLD_HP_Tschebycheff();
int FLD_TP_Tschebycheff();
int FLD_Referenz_TP_Butterworth(void);
int FLD_BS_Butterworth();
int FLD_BP_Butterworth();
int FLD_HP_Butterworth();
int FLD_TP_Butterworth();
void FLD_Init_Filter();
void FLD_Set_Instance(DSP_DATA *dsp_instance);
void FLD_Digitale_Koeffizienten_ermitteln(void);
void FLD_Coefficient_Assignment(void);
int FLD_BS_Transformation(int filt_char);
int FLD_BP_Transformation(int filt_char);
int FLD_Referenz_TP_Cauer();
int FLD_BS_Cauer();
int FLD_BP_Cauer();
int FLD_HP_Cauer();
int FLD_TP_Cauer();
void FLD_Get_Constants();
float FLD_gammaM(float a, float b);
float FLD_u_frq(float omega_s, int nn);
float FLD_u_db(float as, float ad);
float FLD_Sn(float u, float z);
int FLD_BS_Tschebycheff();
int FLD_BP_Tschebycheff();
int FLD_HP_Tschebycheff();
int FLD_TP_Tschebycheff();
int FLD_Referenz_TP_Tschebycheff();
int FLD_BS_Butterworth();
int FLD_BP_Butterworth();
int FLD_HP_Butterworth();
int FLD_TP_Butterworth();
int FLD_Referenz_TP_Butterworth(void);
void FLD_Init_Filter();




#endif