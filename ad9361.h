/*
 * Copyright (C) 2015 Analog Devices, Inc.
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

#ifndef __AD9361_H__
#define __AD9361_H__

#ifdef __cplusplus
extern "C" {
#endif

/* FLAGS */
#define FIXUP_INTERFACE_TIMING	1
#define CHECK_SAMPLE_RATES	2

#ifdef _WIN32
#   ifdef LIBAD9361_EXPORTS
#	define __api __declspec(dllexport)
#   else
#	define __api __declspec(dllimport)
#   endif
#elif __GNUC__ >= 4
#   define __api __attribute__((visibility ("default")))
#else
#   define __api
#endif

struct iio_context;
struct iio_device;

struct filter_design_parameters {
    double Rdata;
    double Fpass;
    double Fstop;
    double caldiv;
    double FIR;
    double HB1;
    double DAC_div;
    const char *Type;
    const char *RxTx;
    double RFbw;
    double converter_rate;
    double PLL_rate;
    double Fcenter;
    double wnom;
    double FIRdBmin;
    double int_FIR;
    double PLL_mult;
    double Apass;
    double Astop;
    double phEQ;
    double HB2;
    double HB3;
    double maxTaps;
};

__api int ad9361_multichip_sync(struct iio_device *master,
		struct iio_device **slaves, unsigned int num_slaves,
		unsigned int flags);

__api int ad9361_fmcomms5_multichip_sync(
		struct iio_context *ctx, unsigned int flags);

__api int ad9361_set_bb_rate(struct iio_device *dev, unsigned long rate);

__api int ad9361_set_trx_fir_enable(struct iio_device *dev, int enable);

__api int ad9361_get_trx_fir_enable(struct iio_device *dev, int *enable);

__api int ad9361_generate_fir_taps(struct filter_design_parameters *parameters,
                                   short *taps, int *num_taps, int *gain);

__api int ad9361_calculate_rf_clock_chain(unsigned long tx_sample_rate,
                                          unsigned long rate_gov,
                                          unsigned long *rx_path_clks,
                                          unsigned long *tx_path_clks);

__api int ad9361_calculate_rf_clock_chain_fdp(struct filter_design_parameters *fdpTX,
                                              struct filter_design_parameters *fdpRX,
                                              unsigned long sample_rate);

__api int ad9361_set_bb_rate_custom_filter_auto(struct iio_device *dev,
                                                unsigned long rate);

__api int ad9361_set_bb_rate_custom_filter_manual(struct iio_device *dev,
                                                  unsigned long rate, unsigned long Fpass,
                                                  unsigned long Fstop, unsigned long wnom_tx,
                                                  unsigned long wnom_rx);


#ifdef __cplusplus
}
#endif

#undef __api

#endif /* __AD9361_H__ */
