# FilterDesign
Digital Filter Design

This project calculate and test filter coefficients. It runs in MPLABX Simulator for PIC32MX795F512. To visulalize the results the DMCI is used. 

1. Declare a Structure my_filter form type DSP_DATA
2. set the pointer to this structure into the filter library: FLD_Set_Instance(&my_filter);
3. set Parameters:

Parameter|value|Comment
-------- | --- | ---------
my_filter.Samplerate | 8000.0 | Nyquist
my_filter.F1_Analog | 400.0 | 1.corner frequency
my_filter.F2_Analog | 800.0 | 2. corner frequency
my_filter.F3_Analog | 1200.0 | 3. corner frequency
my_filter.F4_Analog | 1800.0 | 4. corner frequency
my_filter.as_Blocking_Attenuation | 40.0 | Ripple in Stop Band
my_filter.ad_Passband_Attenuation | 3.0 | Ripple in Pass Band

4. call FLD_Frequency_Transformation(); this will pre-warp the frequencies to compensate the error of the bilinear transformation
5. call FLD_Init_Filter(); to pre initialize the coefficients
6. call the filter design function. For ex. FLD_BP_Butterworth(); which means Bandpass with Butterworth charackteristic
7. call FLD_Coefficient_Assignment(); this transforms the coefficients to the actual biquad coefficients
8. call FLD_InitCoeffsFloat(); this transforms the biquads to the actual used numerator and denumerator coefficients
9. call FLD_InitCoeffsFixpoint(); for same step as in 8. but for fixpoint coefficients

in 

  \FilterDesign\FilterDesign\filter_design.h

  #define INTEGER_PRECISION   27

this Marco defines the fixpoint "resolution". A 16 would mean Q15


10. calculate step response with 


  float FLT_float_sos(float *denumerator, float *numerator, float *taps, float input_value, int sub_filter):

  long FLT_long_sos(long *denumerator, long *numerator, long *taps, long input_value, int sub_filter);


11. calculate FFT in float from float step response with FFT_float_fft()
12. calculate FFT in float from step step response with FFT_float_fft()
13. calculate Absolute values in dB from FFT result for float and fix
14. transfomr the results to fixpoint because the DMCi cannot handle float correctly
15. open DMCI plugin
16. open the dmci_config.xml configuration
17. place a breakpoint at the end of main on the NOP. 
18. run the simulation
19. look in the DMCI windows for the results

![DMCI](/images/dmci.png)

