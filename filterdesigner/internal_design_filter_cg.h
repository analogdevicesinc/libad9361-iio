/*
 * Sponsored Third Party Support License -- for use only to support
 * products interfaced to MathWorks software under terms specified in your
 * company's restricted use license agreement.
 * File: internal_design_filter_cg.h
 *
 * MATLAB Coder version            : 3.2
 * C/C++ source code generated on  : 15-Aug-2017 18:43:56
 */

#ifndef INTERNAL_DESIGN_FILTER_CG_H
#define INTERNAL_DESIGN_FILTER_CG_H

/* Include Files */
#include <float.h>
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "rt_defines.h"
#include "rt_nonfinite.h"
#include "rtwtypes.h"
#include "internal_design_filter_cg_types.h"

/* Type Definitions */
#include <stdio.h>

/* Function Declarations */
#ifdef __cplusplus

extern "C" {

#endif

  extern void internal_design_filter_cg(double Rdata, double Fpass, double Fstop,
    double caldiv, double FIR, double HB1, double PLL_mult, double Apass, double
    Astop, double phEQ, double HB2, double HB3, const char Type[7], const char
    RxTx[2], double RFbw, double DAC_div, double converter_rate, double PLL_rate,
    double Fcenter, double wnom, double FIRdBmin, double int_FIR, short
    outputTaps[128]);
  extern void internal_design_filter_cg_initialize(void);
  extern void internal_design_filter_cg_terminate(void);

#ifdef __cplusplus

}
#endif
#endif

/*
 * File trailer for internal_design_filter_cg.h
 *
 * [EOF]
 */
