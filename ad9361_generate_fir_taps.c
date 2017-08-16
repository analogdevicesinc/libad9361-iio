/*
 * Copyright (C) 2016 Analog Devices, Inc.
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

#include "ad9361.h"
#include <errno.h>
#include <stdio.h>
#include "filterdesigner/internal_design_filter_cg.h"

#ifdef _WIN32
#include <windows.h>
#else
#include <unistd.h>
#endif

#ifdef _MSC_BUILD
#define snprintf sprintf_s
#endif

void ad9361_generate_fir_taps(struct filter_design_parameters *parameters, short *taps)
{
	// Call filter designer
	internal_design_filter_cg_initialize();
	internal_design_filter_cg(parameters->Rdata, parameters->Fpass, parameters->Fstop,
	 	parameters->caldiv, parameters->FIR, parameters->HB1, parameters->PLL_mult,
	 	parameters->Apass, parameters->Astop, parameters->phEQ, parameters->HB2,
	 	parameters->HB3, parameters->Type, parameters->RxTx, parameters->RFbw,
	 	parameters->DAC_div, parameters->converter_rate, parameters->PLL_rate,
	 	parameters->Fcenter, parameters->wnom, parameters->FIRdBmin, parameters->int_FIR,
	 	taps);
	internal_design_filter_cg_terminate();
}
