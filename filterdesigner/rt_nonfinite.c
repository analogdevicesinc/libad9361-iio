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

/*
 * Abstract:
 *      MATLAB for code generation function to initialize non-finites,
 *      (Inf, NaN and -Inf).
 */
#include "rt_nonfinite.h"
#include <math.h>

real_T rtInf;
real_T rtMinusInf;
real_T rtNaN;
real32_T rtInfF;
real32_T rtMinusInfF;
real32_T rtNaNF;

/* Function: rt_InitInfAndNaN ==================================================
 * Abstract:
 * Initialize the rtInf, rtMinusInf, and rtNaN needed by the
 * generated code. NaN is initialized as non-signaling. Assumes IEEE.
 */
/* Suppress Visual Studio 2013 INFINITY macro expansion compiler warning. */
#if defined(_MSC_VER) && _MSC_VER == 1800

#pragma warning(disable: 4756 56)

#endif

void rt_InitInfAndNaN(size_t realSize)
{
  (void)realSize;
  rtNaN = nan("");
  rtNaNF = nanf("");
  rtInf = (real_T)INFINITY;
  rtInfF = (real32_T)INFINITY;
  rtMinusInf = -(real_T)INFINITY;
  rtMinusInfF = -(real32_T)INFINITY;

#if defined(_MSC_VER) && _MSC_VER == 1800

#pragma warning(default: 4756 56)

#endif

}

/* Function: rtIsInf ==================================================
 * Abstract:
 * Test if value is infinite
 */
boolean_T rtIsInf(real_T value)
{
  return (isinf(value) ? 1U : 0U);
}

/* Function: rtIsInfF =================================================
 * Abstract:
 * Test if single-precision value is infinite
 */
boolean_T rtIsInfF(real32_T value)
{
  return (isinf((real_T)value) ? 1U : 0U);
}

/* Function: rtIsNaN ==================================================
 * Abstract:
 * Test if value is not a number
 */
boolean_T rtIsNaN(real_T value)
{
  return (isnan(value) ? 1U : 0U);
}

/* Function: rtIsNaNF =================================================
 * Abstract:
 * Test if single-precision value is not a number
 */
boolean_T rtIsNaNF(real32_T value)
{
  return (isnan((real_T)value) ? 1U : 0U);
}

/*
 * File trailer for rt_nonfinite.c
 *
 * [EOF]
 */
