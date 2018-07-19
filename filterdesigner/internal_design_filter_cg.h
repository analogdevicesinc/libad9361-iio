/*
 * Copyright (C) 2017 Analog Devices, Inc.
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
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

/* Function Declarations */
#ifdef __cplusplus

extern "C" {

#endif

extern void internal_design_filter_cg(double Rdata, double Fpass, double Fstop,
                                      double caldiv, double FIR, double HB1, double PLL_mult, double Apass, double
                                      Astop, double phEQ, double HB2, double HB3, const char Type[7], const char
                                      RxTx[2], double RFbw, double DAC_div, double converter_rate, double PLL_rate,
                                      double Fcenter, double wnom, double FIRdBmin, double int_FIR, double maxTaps,
                                      short outputTaps[128], double *numOutputTaps, double *filterGain);
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
