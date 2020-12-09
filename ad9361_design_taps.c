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
#include "filterdesigner/internal_design_filter_cg.h"
#include <errno.h>
#include <iio.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>

#ifdef _WIN32
#include <windows.h>
#else
#include <unistd.h>
#endif
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#define FIR_BUF_SIZE	8192

#ifdef _MSC_BUILD
#define snprintf sprintf_s
#endif

int ad9361_generate_fir_taps(struct filter_design_parameters *parameters,
                             short *taps, int *num_taps, int *gain)
{
    double dnum_taps = 0;
    double dgain = 0;

    // Call filter designer
    internal_design_filter_cg_initialize();
    // Designer will alway return a filter, but it may not meet specifications
    internal_design_filter_cg(
        parameters->Rdata, parameters->Fpass, parameters->Fstop,
        parameters->caldiv, parameters->FIR, parameters->HB1,
        parameters->PLL_mult, parameters->Apass, parameters->Astop,
        parameters->phEQ, parameters->HB2, parameters->HB3, parameters->Type,
        parameters->RxTx, parameters->RFbw, parameters->DAC_div,
        parameters->converter_rate, parameters->PLL_rate, parameters->Fcenter,
        parameters->wnom, parameters->FIRdBmin, parameters->int_FIR,
        parameters->maxTaps, taps, &dnum_taps, &dgain);
    internal_design_filter_cg_terminate();
    *num_taps = (int)dnum_taps;
    *gain = (int)dgain;
    // Filter with less than 32 taps is an error
    if (*num_taps < 32)
        return -EDOM;
    else
        return 0;
}


double calculate_rfbw(double pll_rate, double caldiv, bool TX,
                      double *rcaldiv)
{
    double rfbw, min_rfbw, max_rfbw, scale;
    if (TX) {
        scale = 1.6;
        min_rfbw = 1250000;
        max_rfbw = 40000000;
    } else {
        scale = 1.4;
        min_rfbw = 400000;
        max_rfbw = 56000000;
    }
    rfbw =
        (double)round((pll_rate / caldiv) * (2 / (scale * (2 * M_PI) / log(2))));

    // If the RF bandwidth is outside the range of acceptable values we modify
    // the divider value until it falls into an acceptable range.
    while ((rfbw < min_rfbw) || (rfbw > max_rfbw)) {
        if (rfbw < min_rfbw)
            caldiv = caldiv - 1;
        else
            caldiv = caldiv + 1;

        if ((caldiv < 1) || (caldiv > 511)) {
            fprintf(stderr,"Calibration divider out of bounds (1 - 511): %f\n", caldiv);
            return -EINVAL;
        }
        rfbw = calculate_rfbw(pll_rate, caldiv, TX, rcaldiv);
    }
    *rcaldiv = caldiv;
    return rfbw;
}

void set_max_taps(struct filter_design_parameters *fdpTX,
                  struct filter_design_parameters *fdpRX)
{
    // RX side
    int N,M,K;
    if (fdpRX->HB3 == 3)
        N = 16*floor(fdpRX->converter_rate/(fdpRX->Rdata));
    else
        N = 16*floor(fdpRX->converter_rate/(2*fdpRX->Rdata));
    if (N>128)
        N = 128;
    // TX side
    if (fdpTX->FIR==1)
        M = 64;
    else
        M = 128;
    K = 16*floor(fdpTX->converter_rate*fdpTX->DAC_div/(2*fdpTX->Rdata));
    if (K<M)
        M = K;

    // Pick the smallest
    if (M>N) {
        fdpTX->maxTaps = N;
        fdpRX->maxTaps = N;
    } else {
        fdpTX->maxTaps = M;
        fdpRX->maxTaps = M;
    }
}

int build_configuration(struct filter_design_parameters *fdpTX,
                        struct filter_design_parameters *fdpRX,
                        unsigned long sample_rate,
                        unsigned long Fpass,
                        unsigned long Fstop,
                        unsigned long wnomTX,
                        unsigned long wnomRX)
{

    double div, max;
    unsigned long rx_path_clk[6];
    unsigned long tx_path_clk[6];
    unsigned long *path_clk;
    struct filter_design_parameters *fdp;
    int ret,k;
    unsigned long rate_gov = 0;

    ret = ad9361_calculate_rf_clock_chain((unsigned long) sample_rate,
                                          rate_gov, rx_path_clk, tx_path_clk);
    if (ret < 0)
        return -EINVAL;

    for (k=0; k<2; k++) {

        if (k > 0) {
            path_clk = tx_path_clk;
            fdp = fdpTX;
            fdp->RxTx = "Tx";
            fdp->DAC_div = (double) rx_path_clk[1]/tx_path_clk[1];
        } else {
            path_clk = rx_path_clk;
            fdp = fdpRX;
            fdp->RxTx = "Rx";
            fdp->DAC_div = 1.0;
        }
        // Map rates and dividers
        fdp->PLL_rate = (double) path_clk[0];
        fdp->converter_rate = (double) path_clk[1];
        fdp->PLL_mult = (double) path_clk[0]/path_clk[1];
        fdp->HB3 = (double) path_clk[1]/path_clk[2];
        fdp->HB2 = (double) path_clk[2]/path_clk[3];
        fdp->HB1 = (double) path_clk[3]/path_clk[4];
        fdp->FIR = (double) path_clk[4]/path_clk[5];

        // Set default parameters
        fdp->Rdata = (double)path_clk[5];
        fdp->Type = "Lowpass";
        fdp->int_FIR = 1;
        fdp->Apass = 0.5;
        fdp->Astop = 80;
        fdp->phEQ = -1;
        fdp->FIRdBmin = 0;
        // Define filter design specifications
        fdp->Fpass = (double) Fpass;
        fdp->Fstop = (double) Fstop;
        fdp->Fcenter = 0.0;
        if (k>0)
            fdp->wnom = (double) wnomTX;
        else
            fdp->wnom = (double) wnomRX;
        // Determine default analog bandwidth
        div = ceil((fdp->PLL_rate / fdp->wnom) * (log(2) / (2 * M_PI)));
        max = (div > 1) ? div : 1.0;
        fdp->caldiv = (max > 511) ? 511.0 : max;
        fdp->RFbw = calculate_rfbw(fdp->PLL_rate,fdp->caldiv, k>0, &(fdp->caldiv));

        if (fdp->RFbw < 0)
            return -EINVAL;
    }
    set_max_taps(fdpTX,fdpRX);

    return 0;

}

int apply_custom_filter(struct iio_device *dev, unsigned dec_tx,
                        unsigned dec_rx, short *tapsTx,
                        short *tapsRx, unsigned taps,
                        unsigned long rate,
                        int gain_tx, int gain_rx,
                        unsigned long wnom_tx, unsigned long wnom_rx)
{
    struct iio_channel *chanTX, *chanRX;
    long long current_rate;
    int ret, i, enable, len = 0;
    char *buf;

    chanTX = iio_device_find_channel(dev, "voltage0", true);
    if (chanTX == NULL)
        return -ENODEV;

    ret = iio_channel_attr_read_longlong(chanTX, "sampling_frequency", &current_rate);
    if (ret < 0)
        return ret;

    ret = ad9361_get_trx_fir_enable(dev, &enable);
    if (ret < 0)
        return ret;

    if (enable) {
        if (current_rate <= (25000000 / 12))
            iio_channel_attr_write_longlong(chanTX, "sampling_frequency", 3000000);

        ret = ad9361_set_trx_fir_enable(dev, false);
        if (ret < 0)
            return ret;
    }

    buf = (char*) malloc(FIR_BUF_SIZE);
    if (!buf)
        return -ENOMEM;


    len += snprintf(buf + len, FIR_BUF_SIZE - len, "RX 3 GAIN %d DEC %d\n", gain_rx,
                    dec_rx);
    len += snprintf(buf + len, FIR_BUF_SIZE - len, "TX 3 GAIN %d INT %d\n", gain_tx,
                    dec_tx);

    for (i = 0; i < (int)taps; i++)
        len += snprintf(buf + len, FIR_BUF_SIZE - len, "%d,%d\n", tapsRx[i], tapsTx[i]);

    len += snprintf(buf + len, FIR_BUF_SIZE - len, "\n");

    ret = iio_device_attr_write_raw(dev, "filter_fir_config", buf, len);
    free (buf);

    if (ret < 0)
        return ret;

    if (rate <= (25000000 / 12))  {

        int dacrate, txrate, max;
        char readbuf[100];
        ret = iio_device_attr_read(dev, "tx_path_rates", readbuf, sizeof(readbuf));
        if (ret < 0)
            return ret;
        ret = sscanf(readbuf, "BBPLL:%*d DAC:%d T2:%*d T1:%*d TF:%*d TXSAMP:%d",
                     &dacrate, &txrate);
        if (ret != 2)
            return -EFAULT;
        if (txrate == 0)
            return -EINVAL;
        max = (dacrate / txrate) * 16;
        if (max < taps)
            iio_channel_attr_write_longlong(chanTX, "sampling_frequency", 3000000);


        ret = ad9361_set_trx_fir_enable(dev, true);
        if (ret < 0)
            return ret;
        ret = iio_channel_attr_write_longlong(chanTX, "sampling_frequency", rate);
        if (ret < 0)
            return ret;
    } else {
        ret = iio_channel_attr_write_longlong(chanTX, "sampling_frequency", rate);
        if (ret < 0)
            return ret;
        ret = ad9361_set_trx_fir_enable(dev, true);
        if (ret < 0)
            return ret;
    }

    chanRX = iio_device_find_channel(dev, "voltage0", false);
    if (chanRX == NULL)
        return -ENODEV;
    ret = iio_channel_attr_write_longlong(chanTX, "rf_bandwidth", wnom_tx);
    if (ret < 0)
        return ret;
    ret = iio_channel_attr_write_longlong(chanRX, "rf_bandwidth", wnom_rx);
    if (ret < 0)
        return ret;

    return 0;
}

int ad9361_calculate_rf_clock_chain_fdp(struct filter_design_parameters *fdpTX,
                                        struct filter_design_parameters *fdpRX,
                                        unsigned long sample_rate)
{
    int ret;

    // Set default configuration
    unsigned long Fpass = sample_rate / 3.0;
    unsigned long Fstop = Fpass * 1.25;
    unsigned long wnomTX = 1.6 * Fstop;
    unsigned long wnomRX = 1.4 * Fstop;

    ret = build_configuration(fdpTX, fdpRX, sample_rate, Fpass, Fstop, wnomTX,
                              wnomRX);
    if (ret<0)
        return ret;

    return 0;

}

int ad9361_set_bb_rate_custom_filter_auto(struct iio_device *dev,
        unsigned long rate)
{
    struct filter_design_parameters fdpTX;
    struct filter_design_parameters fdpRX;
    short taps_tx[128];
    short taps_rx[128];
    int ret, num_taps_tx, num_taps_rx, gain_tx, gain_rx;
    unsigned dec_tx, dec_rx, num_taps;

    ret = ad9361_calculate_rf_clock_chain_fdp(&fdpTX, &fdpRX, rate);
    if (ret < 0)
        return ret;

    ret = ad9361_generate_fir_taps(&fdpRX, taps_rx, &num_taps_rx, &gain_rx);
    if (ret < 0)
        return ret;

    ret = ad9361_generate_fir_taps(&fdpTX, taps_tx, &num_taps_tx, &gain_tx);
    if (ret < 0)
        return ret;

    dec_tx = (unsigned) fdpTX.FIR;
    dec_rx = (unsigned) fdpRX.FIR;
    num_taps = (unsigned) fdpTX.maxTaps;
    ret = apply_custom_filter(dev, dec_tx, dec_rx, taps_tx, taps_rx, num_taps,
                              rate, gain_tx, gain_rx, fdpTX.wnom, fdpRX.wnom);
    if (ret < 0)
        return ret;

    return 0;
}

int ad9361_set_bb_rate_custom_filter_manual(struct iio_device *dev,
        unsigned long rate, unsigned long Fpass,
        unsigned long Fstop, unsigned long wnom_tx, unsigned long wnom_rx)
{
    struct filter_design_parameters fdpTX;
    struct filter_design_parameters fdpRX;
    short taps_tx[128];
    short taps_rx[128];
    int ret, num_taps_tx, num_taps_rx, gain_tx, gain_rx;
    unsigned dec_tx, dec_rx, num_taps;

    if (Fpass >= Fstop)
        return -EINVAL;

    ret = build_configuration(&fdpTX, &fdpRX, rate, Fpass, Fstop, wnom_tx,
                              wnom_rx);
    if (ret<0)
        return ret;

    ret = ad9361_generate_fir_taps(&fdpRX, taps_rx, &num_taps_rx, &gain_rx);
    if (ret < 0)
        return ret;

    ret = ad9361_generate_fir_taps(&fdpTX, taps_tx, &num_taps_tx, &gain_tx);
    if (ret < 0)
        return ret;

    dec_tx = (unsigned) fdpTX.FIR;
    dec_rx = (unsigned) fdpRX.FIR;
    num_taps = (unsigned) fdpTX.maxTaps;

    ret = apply_custom_filter(dev, dec_tx, dec_rx, taps_tx, taps_rx, num_taps,
                              rate, gain_tx, gain_rx, wnom_tx, wnom_rx);
    if (ret < 0)
        return ret;

    return 0;
}
