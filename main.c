
#include <math.h>
#include <float.h>
#include <limits.h>
#include <stdint.h>

#include ".\alg\fft_float.h"
#include ".\alg\filter_sos.h"
#include ".\FilterDesign\filter_design.h"
#include "alg/fft_long.h"

#pragma config CP = OFF, BWP = OFF, PWP = OFF
#pragma config FNOSC = FRCPLL, POSCMOD = OFF
#pragma config FPLLIDIV = DIV_2, FPLLMUL = MUL_18, FPBDIV = DIV_2, FPLLODIV = DIV_1
#pragma config FWDTEN = OFF, FSOSCEN = OFF, IESO = OFF, OSCIOFNC = OFF

#define NOP asm("nop")

long long_round_float(float x);

#define MAX_BIQADS	16
#define MAX_BQ_SIZE (MAX_BIQADS*3) 

#define SIZE_OF_DATA	512
#define SIZE_OF_RESULT	(SIZE_OF_DATA>>1)

float_complex float_complex_FloatData[SIZE_OF_DATA];
long_complex long_complex_FloatData[SIZE_OF_DATA];

float_complex float_complex_FixpointData[SIZE_OF_DATA];
long_complex long_complex_FixpointData[SIZE_OF_DATA];

float float_StepResponse_FloatData[SIZE_OF_DATA];
float float_FFTResult_FloatData[SIZE_OF_RESULT];
float float_FFTResult_FixpointData[SIZE_OF_RESULT];
float float_FFTResult_FFTFixpointData[SIZE_OF_RESULT];

long long_StepResponse_FloatData[SIZE_OF_DATA]; // DMCI 1/1
long long_FFTResult_FloatData[SIZE_OF_RESULT]; // DMCI 2/1
long long_StepResponse_FixpointData[SIZE_OF_DATA]; // DMCI 1/2
long long_FFTResult_FixpointData[SIZE_OF_RESULT]; // DMCI 2/2

float float_Numerator[MAX_BQ_SIZE];
float float_Denumerator[MAX_BQ_SIZE];
float float_Taps[MAX_BQ_SIZE];

long long_Numerator[MAX_BQ_SIZE];
long long_Denumenrator[MAX_BQ_SIZE];
long long_Taps[MAX_BQ_SIZE];

DSP_DATA my_filter;

int main(void) {
    int i;
    float max_val;

    //--------------------------------------------------------------------------
    FLD_Set_Instance(&my_filter);
    //--------------------------------------------------------------------------
    my_filter.Samplerate = 8000.0;
    my_filter.F1_Analog = 400.0; // 1. corner frequency
    my_filter.F2_Analog = 800.0; // 2. corner frequency
    my_filter.F3_Analog = 1200.0; // 3. corner frequency
    my_filter.F4_Analog = 1800.0; // 4. corner frequency
    my_filter.as_Blocking_Attenuation = 40.0; // Ripple in Stop Band
    my_filter.ad_Passband_Attenuation = 3.0; // Ripple in Pass Band
    //--------------------------------------------------------------------------
    FLD_Frequency_Transformation();
    FLD_Init_Filter();
    //--------------------------------------------------------------------------
    //FLD_TP_Butterworth();
    //FLD_HP_Butterworth();
    FLD_BP_Butterworth();
    //FLD_BS_Butterworth();
    //--------------------------------------------------------------------------
    //FLD_TP_Tschebycheff();
    //FLD_HP_Tschebycheff();
    //FLD_BP_Tschebycheff();
    //FLD_BS_Tschebycheff();
    //--------------------------------------------------------------------------
    //FLD_TP_Cauer();
    //FLD_HP_Cauer();
    FLD_BP_Cauer();
    //FLD_BS_Cauer();
    //--------------------------------------------------------------------------
    FLD_Coefficient_Assignment();
    FLD_InitCoeffsFloat(float_Numerator, float_Denumerator, float_Taps, my_filter.Sub_Filter);
    FLD_InitCoeffsFixpoint(long_Numerator, long_Denumenrator, long_Taps, my_filter.Sub_Filter);
    //--------------------------------------------------------------------------       
    // Data Array Zeroing
    for (i = 0; i < SIZE_OF_DATA; i++) {
        float_complex_FloatData[i].x = 0.0;
        float_complex_FloatData[i].y = 0.0;
        long_complex_FloatData[i].x = 0;
        long_complex_FloatData[i].y = 0;
        float_complex_FixpointData[i].x = 0.0;
        float_complex_FixpointData[i].y = 0.0;
        long_complex_FixpointData[i].x = 0;
        long_complex_FixpointData[i].y = 0;
    }
    //==========================================================================



    //==========================================================================
    // Filter calculate step response in ==>> Float <<==
    float_complex_FloatData[0].x = 1.0;
    for (i = 0; i < SIZE_OF_DATA; i++) {
        float_StepResponse_FloatData[i] = float_complex_FloatData[i].x =
                FLT_float_sos(float_Denumerator, float_Numerator, float_Taps,
                float_complex_FloatData[i].x, my_filter.Sub_Filter);

        long_StepResponse_FloatData[i] =
                long_round_float(float_StepResponse_FloatData[i] * (float) QX_ONE);

        long_complex_FloatData[i].x = long_StepResponse_FloatData[i];
    }
    //--------------------------------------------------------------------------                  
    // FFT of Float calculation
    FFT_float_fft(SIZE_OF_DATA, float_complex_FloatData);
    //--------------------------------------------------------------------------
    // FFT Data to dB
    // This is necessary because the DMCI cannot display float correctly
    for (i = 0; i < SIZE_OF_RESULT; i++) {
        float_FFTResult_FloatData[i] =
                20 * log10f(sqrt(pow(float_complex_FloatData[i].x, 2)
                + pow(float_complex_FloatData[i].y, 2)));

        long_FFTResult_FloatData[i] = long_round_float(float_FFTResult_FloatData[i]);
    }
    long_FFTResult_FloatData[0] = long_FFTResult_FloatData[1];
    //==========================================================================





    //==========================================================================
    // Filter caculate step response in ==>> Fixpoint <<==
    long_StepResponse_FixpointData[0] = QX_ONE;
    for (i = 0; i < SIZE_OF_DATA; i++) {
        long_StepResponse_FixpointData[i] =
                FLT_long_sos(long_Denumenrator, long_Numerator, long_Taps,
                long_StepResponse_FixpointData[i],
                my_filter.Sub_Filter);
        float_complex_FixpointData[i].x = (float) long_StepResponse_FixpointData[i];
    }
    //--------------------------------------------------------------------------
    // FFT of Fixpoint calculation
    FFT_float_fft(SIZE_OF_DATA, float_complex_FixpointData);
    //--------------------------------------------------------------------------
    // FFT Data to dB
    // This is necessary because the DMCI cannot display float correctly
    for (i = 0; i < SIZE_OF_RESULT; i++) {
        float_FFTResult_FloatData[i] =
                20 * log10f(sqrt(pow(float_complex_FixpointData[i].x, 2)
                + pow(float_complex_FixpointData[i].y, 2)));

        long_FFTResult_FixpointData[i] = long_round_float(float_FFTResult_FloatData[i]);
    }
    // avoid first data as unvalid becuase of log(0)
    long_FFTResult_FixpointData[0] = long_FFTResult_FixpointData[1];
    // scale to zero
    max_val = -10000.0;
    for (i = 0; i < SIZE_OF_RESULT; i++) {
        if (max_val < long_FFTResult_FixpointData[i]) {
            max_val = long_FFTResult_FixpointData[i];
        }
    }
    for (i = 0; i < SIZE_OF_RESULT; i++) {
        long_FFTResult_FixpointData[i] -= max_val;
    }
    //==========================================================================

    NOP;
}

long long_round_float(float x) {
    float res;

    if (x >= 0.0F) {
        res = ceilf(x);
        if (res - x > 0.5F)
            res -= 1.0F;
    } else {
        res = ceilf(-x);
        if (res + x > 0.5F)
            res -= 1.0F;
        res = -res;
    }
    return (long) res;
}
