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

#include "ad9361.h"
#include <math.h>
#include <stdbool.h>
#include <errno.h>
#include <stdio.h>
#include "filterdesigner/internal_design_filter_cg.h"

#ifdef _WIN32
#include <windows.h>
#define M_PI 3.14159265358979323846
#else
#include <unistd.h>
#endif

#ifdef _MSC_BUILD
#define snprintf sprintf_s
#endif

#define AD9361_MAX_RATE 61440000
#define AD9361_MIN_RATE 520833
#define AD9361_MAX_RX_HB1  245760000
#define AD9361_MAX_RX_HB2  320000000
#define AD9361_MAX_RX_HB3  640000000
#define AD9361_MAX_TX_HB1  160000000
#define AD9361_MAX_TX_HB2  320000000
#define AD9361_MAX_TX_HB3  320000000
#define AD9361_MAX_FIR AD9361_MAX_RATE*2
#define AD9361_MAX_BBPLL_FREQ 1430000000
#define AD9361_MIN_BBPLL_FREQ 715000000
#define AD9361_MIN_ADC_CLK 25000000
#define AD9361_MAX_ADC_CLK 640000000
#define AD9361_MIN_DAC_CLK 25000000
#define AD9361_MAX_DAC_CLK AD9361_MAX_ADC_CLK/2

void ad9361_generate_fir_taps(struct filter_design_parameters *parameters, short *taps, int *num_taps)
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
	// Set number of taps used
	*num_taps = 0;
	for((*num_taps)=0;((*num_taps)<128)&&taps[*num_taps]!=0;(*num_taps)++);
}

int determine_HB_rate(unsigned *options, unsigned num_options,
  unsigned long max_rate, unsigned long min_rate, unsigned long multiplier)
{
  int o=1;
  for (;o<(num_options-1);o++)
  {
    if ((max_rate>=(multiplier*options[num_options-o]))
     && (min_rate<=(multiplier*options[num_options-o])))
      break;
  }
  return options[num_options-o];
}

void define_rates(struct filter_design_parameters *fdp, bool TX, bool dec3enable)
{
  unsigned long HB3, HB2, HB1;
  if (TX)
  {
    fdp->DAC_div = 2;
    HB3 = AD9361_MAX_TX_HB3;
    HB2 = AD9361_MAX_TX_HB2;
    HB1 = AD9361_MAX_TX_HB1;
  }
  else
  {
    fdp->DAC_div = 1;
    HB3 = AD9361_MAX_RX_HB3;
    HB2 = AD9361_MAX_RX_HB2;
    HB1 = AD9361_MAX_RX_HB1;
  }
  // HB3
  unsigned options3[3] = {1, 2, 3};
  unsigned options4[3] = {1, 2, 4};
  unsigned options64[7] = {1, 2, 4, 8, 16, 32, 64};
  if (dec3enable)
    fdp->HB3 = (double) determine_HB_rate(options3,3,HB3,0,fdp->Rdata);
  else
    fdp->HB3 = (double) determine_HB_rate(options3,2,HB3,0,fdp->Rdata);
  fdp->HB2 = (double) determine_HB_rate(options4,2,HB2,0,fdp->Rdata*fdp->HB3);
  fdp->HB1 = (double) determine_HB_rate(options4,2,HB1,0,fdp->Rdata*fdp->HB3*fdp->HB2);
  fdp->FIR = (double) determine_HB_rate(options4,3,AD9361_MAX_FIR,0,fdp->Rdata*fdp->HB3*fdp->HB2*fdp->HB1);
  fdp->PLL_mult = (double) determine_HB_rate(options64,7,AD9361_MAX_BBPLL_FREQ,
    AD9361_MIN_BBPLL_FREQ,fdp->Rdata*fdp->HB3*fdp->HB2*fdp->HB1*fdp->FIR*fdp->DAC_div);
  fdp->converter_rate = fdp->Rdata * fdp->HB3 * fdp->HB2 * fdp->HB1 * fdp->FIR;
  fdp->PLL_rate = fdp->converter_rate * fdp->DAC_div * fdp->PLL_mult;
}

double calculate_rfbw(double pll_rate, double caldiv, bool TX, double *rcaldiv)
{
  double rfbw, min_rfbw, max_rfbw, *tmp, scale;
  if (TX)
  {
    scale = 1.6;
    min_rfbw = 1250000;
    max_rfbw = 40000000;
  }
  else
  {
    scale = 1.4;
    min_rfbw = 400000;
    max_rfbw = 56000000;
  }
  rfbw = (double) round((pll_rate/ caldiv)*(2/(scale*(2*M_PI)/log(2))));

  // If the RF bandwidth is outside the range of acceptable values we modify
  // the divider value until it falls into an acceptable range.
  while ((rfbw < min_rfbw) || (rfbw > max_rfbw))
  {
    if (rfbw < min_rfbw)
        caldiv = caldiv - 1;
    else
        caldiv = caldiv + 1;

    if ((caldiv < 1) || (caldiv > 511))
    {
        printf("Calibration divider out of bounds (1 - 511): %f\n", caldiv);
        return -1;
    }
    rfbw = calculate_rfbw(pll_rate, caldiv, TX, rcaldiv);
  }
  *rcaldiv = caldiv;
  return rfbw;
}

int ad9361_filter_config_from_rate(struct filter_design_parameters *fdp,	unsigned long rate, int TX)
{
  // Check input rate is valid
  if ((rate>AD9361_MAX_RATE) || (rate<AD9361_MIN_RATE))
    return -1;
  // Set default parameters
  fdp->Rdata = (double) rate;
  fdp->Type = "Lowpass";
  if (TX>0)
    fdp->RxTx = "Tx";
  else
    fdp->RxTx = "Rx";
  fdp->int_FIR = 1;
  fdp->Apass = 0.5;
  fdp->Astop = 80;
  fdp->phEQ = -1;
  fdp->FIRdBmin = 0;
	// Determine path rates
  define_rates(fdp,TX,false);
  // Check if rate bounds are met
  bool pll_oob = (fdp->PLL_rate>AD9361_MAX_BBPLL_FREQ) || (fdp->PLL_rate<AD9361_MIN_BBPLL_FREQ);
  bool converter_oob;
  if (TX)
    converter_oob = (fdp->converter_rate > AD9361_MAX_DAC_CLK) || (fdp->converter_rate < AD9361_MIN_DAC_CLK);
  else
    converter_oob = (fdp->converter_rate > AD9361_MAX_ADC_CLK) || (fdp->converter_rate < AD9361_MIN_ADC_CLK);
  // Run with additional decimation option in HB3
  if (converter_oob || pll_oob)
    define_rates(fdp,TX,true);
	// Define filter design specifications
  fdp->Fpass = fdp->Rdata / 3.0;
  fdp->Fstop = fdp->Fpass * 1.25;
  fdp->Fcenter = 0.0;
  if (TX)
    fdp->wnom = 1.6 * fdp->Fstop;
  else
    fdp->wnom = 1.4 * fdp->Fstop;
	// Determine analog bandwidth
  double div = ceil((fdp->PLL_rate/ fdp->wnom)*(log(2)/(2*M_PI)));
  double max = (div>1) ? div : 1.0;
  fdp->caldiv = (max>511) ? 511.0 : max;
  fdp->RFbw = calculate_rfbw(fdp->PLL_rate, fdp->caldiv, TX>0, &(fdp->caldiv));

  return -((int)(fdp->RFbw<0));
}
