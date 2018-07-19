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
 *       MATLAB for code generation function to initialize non-finite, NaN
 */
#include "rtGetNaN.h"
#define NumBitsPerChar                 8U

/* Function: rtGetNaN ==================================================
 * Abstract:
 * Initialize rtNaN needed by the generated code.
 * NaN is initialized as non-signaling. Assumes IEEE.
 */
real_T rtGetNaN(void)
{
    real_T nan = 0.0;
    uint16_T one = 1U;
    enum {
        LittleEndian,
        BigEndian
    } machByteOrder = (*((uint8_T *) &one) == 1U) ? LittleEndian : BigEndian;
    switch (machByteOrder) {
    case LittleEndian: {
        union {
            LittleEndianIEEEDouble bitVal;
            real_T fltVal;
        } tmpVal;

        tmpVal.bitVal.words.wordH = 0xFFF80000U;
        tmpVal.bitVal.words.wordL = 0x00000000U;
        nan = tmpVal.fltVal;
        break;
    }

    case BigEndian: {
        union {
            BigEndianIEEEDouble bitVal;
            real_T fltVal;
        } tmpVal;

        tmpVal.bitVal.words.wordH = 0x7FFFFFFFU;
        tmpVal.bitVal.words.wordL = 0xFFFFFFFFU;
        nan = tmpVal.fltVal;
        break;
    }
    }

    return nan;
}

/* Function: rtGetNaNF ==================================================
 * Abstract:
 * Initialize rtNaNF needed by the generated code.
 * NaN is initialized as non-signaling. Assumes IEEE.
 */
real32_T rtGetNaNF(void)
{
    IEEESingle nanF = { { 0 } };

    uint16_T one = 1U;
    enum {
        LittleEndian,
        BigEndian
    } machByteOrder = (*((uint8_T *) &one) == 1U) ? LittleEndian : BigEndian;
    switch (machByteOrder) {
    case LittleEndian: {
        nanF.wordL.wordLuint = 0xFFC00000U;
        break;
    }

    case BigEndian: {
        nanF.wordL.wordLuint = 0x7FFFFFFFU;
        break;
    }
    }

    return nanF.wordL.wordLreal;
}

/*
 * File trailer for rtGetNaN.c
 *
 * [EOF]
 */
