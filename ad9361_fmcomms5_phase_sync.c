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
#include <iio.h>
#include <stdio.h>
#ifdef _WIN32
#include <windows.h>
#else
#include <unistd.h>
#endif

#include <stdbool.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

#ifdef _MSC_BUILD
#define snprintf sprintf_s
#endif

// Device names
#define DEV_RX_NAME "cf-ad9361-A"
#define DEV_RX_SLAVE_NAME "cf-ad9361-B"
#define DEV_TX_NAME "cf-ad9361-dds-core-lpc"
#define DEV_TX_SLAVE_NAME "cf-ad9361-dds-core-B"
#define DEV_PHY_NAME "ad9361-phy"
#define DEV_PHY_SLAVE_NAME "ad9361-phy-B"

#define DDS_SCALE 0.3
#define SAMPLES 16384
#define TOLERANCE 0.001
#define CALIBRATE_TRIES 100
#define STEP_SIZE 0.1
#define M_2PI 2*M_PI
#define DEBUG 0

#define CHECK(expr) if (expr<0) {return expr;}

static struct iio_device *dev_phy, *dev_phy_slave;
static struct iio_device *dev_rx, *dev_rx_slave;
static struct iio_device *dev_tx, *dev_tx_slave;
static struct iio_channel *dds_out[2][8];
static struct iio_buffer  *rxbuf;
static struct iio_channel *rxa_chan_real, *rxa_chan_imag;
static struct iio_channel *rxb_chan_real, *rxb_chan_imag;

double scale_phase_0_360(double val)
{
	if (val >= 360.0)
		val -= 360.0;

	if (val < 0)
		val += 360.0;

	return val;
}

void dds_tx_phase_rotation(struct iio_device *dev, double val)
{
	long long i, q;
	int d, j;

	if (dev == dev_tx_slave)
		d = 1;
	else
		d = 0;

	i = scale_phase_0_360(val + 90.0) * 1000;
	q = scale_phase_0_360(val) * 1000;

	for (j = 0; j < 8; j++) {
		switch (j) {
			case 0:
			case 1:
			case 4:
			case 5:
				iio_channel_attr_write_longlong(dds_out[d][j], "phase", i);
				break;
			default:
				iio_channel_attr_write_longlong(dds_out[d][j], "phase", q);
		}
	}
}

double unwrap(double x)
{
    x = fmod(x + M_PI, M_2PI);
    if (x < 0)
        x += M_2PI;
    return x - M_PI;
}

double calculate_phase(int16_t *chan0_r, int16_t *chan0_i, int16_t *chan1_r, int16_t *chan1_i, int samples)
{
  double angle_wrapped, angle = 0, prev = 0;
  int k = 0;
  for (;k<samples;k++)
  {
     angle_wrapped = atan2((double)chan1_i[k],(double)chan1_r[k])
                   - atan2((double)chan0_i[k],(double)chan0_r[k]);
    angle += (prev + unwrap(angle_wrapped - prev));
  }
  return angle/k;
}

void near_end_loopback_ctrl(unsigned channel, bool enable)
{
	unsigned tmp;
	struct iio_device *dev = (channel > 3) ?
		dev_rx : dev_rx_slave;
	if (!dev)
		return;

	if (channel > 3)
		channel -= 4;

	if (iio_device_reg_read(dev, 0x80000418 + channel * 0x40, &tmp))
		return;

	if (enable)
		tmp |= 0x1;
	else
		tmp &= ~0xF;

	iio_device_reg_write(dev, 0x80000418 + channel * 0x40, tmp);
}

void configure_ports(unsigned val)
{
	unsigned lp_slave, lp_master, sw;
	char *rx_port, *tx_port;

  // https://wiki.analog.com/resources/eval/user-guides/ad-fmcomms5-ebz/multi-chip-sync#rf_phase_difference

	/*
	*  0 DISABLE: Use RF ports
	*  1 TX1B_B (HPC) -> RX1C_B (HPC) : BIST_LOOPBACK on A
	*  2 TX1B_A (LPC) -> RX1C_B (HPC) : BIST_LOOPBACK on A
	*  3 TX1B_B (HPC) -> RX1C_A (LPC) : BIST_LOOPBACK on B
	*  4 TX1B_A (LPC) -> RX1C_A (LPC) : BIST_LOOPBACK on B
	*
	*/
	switch (val) {
	default:
	case 0:
		lp_slave = 0;
		lp_master = 0;
		sw = 0;
		tx_port = "A";
		rx_port = "A_BALANCED";
		break;
	case 1:
	case 2://
		lp_slave = 0;
		lp_master = 1;
		sw = val - 1;
		tx_port = "B";
		rx_port = "C_BALANCED";
		break;
	case 3:
	case 4:
		lp_slave = 1;
		lp_master = 0;
		sw = val - 1;
		tx_port = "B";
		rx_port = "C_BALANCED";
		break;
	}

	// Set up ports for FPGA BIST Loopback
	near_end_loopback_ctrl(0, lp_slave); /* HPC */
	near_end_loopback_ctrl(1, lp_slave); /* HPC */
	near_end_loopback_ctrl(4, lp_master); /* LPC */
	near_end_loopback_ctrl(5, lp_master); /* LPC */

  // Configure ADG918 switches
	iio_device_debug_attr_write_longlong(dev_phy, "calibration_switch_control", sw);

  // Map ports to switch orientation
	iio_channel_attr_write(iio_device_find_channel(dev_phy, "voltage0", false),"rf_port_select", rx_port);
	iio_channel_attr_write(iio_device_find_channel(dev_phy, "voltage0", true),"rf_port_select", tx_port);
	iio_channel_attr_write(iio_device_find_channel(dev_phy_slave, "voltage0", false),"rf_port_select", rx_port);
	iio_channel_attr_write(iio_device_find_channel(dev_phy_slave, "voltage0", true),"rf_port_select", tx_port);
}

int trx_phase_rotation(struct iio_device *dev, double val)
{
	struct iio_channel *out0, *out1;
	double phase, vcos, vsin;
	unsigned offset;

	bool output = (dev == dev_tx_slave) || (dev == dev_tx);

	phase = val * 2 * M_PI / 360.0;

	vcos = cos(phase);
	vsin = sin(phase);

	if (output)  {
		double corr;
		corr = 1.0 / fmax(fabs(sin(phase) + cos(phase)),
				  fabs(cos(phase) - sin(phase)));
		vcos *= corr;
		vsin *= corr;
	}

	/* Set both RX1 and RX2 */
	for (offset = 0; offset <= 2; offset += 2) {
		if (offset == 2) {
			out0 = iio_device_find_channel(dev, "voltage2", output);
			out1 = iio_device_find_channel(dev, "voltage3", output);
		} else {
			out0 = iio_device_find_channel(dev, "voltage0", output);
			out1 = iio_device_find_channel(dev, "voltage1", output);
		}
		if ((out0 == NULL) || (out0 == NULL))
			return -ENODEV;

		if (out1 && out0) {
			iio_channel_attr_write_double(out0, "calibscale", (double) vcos);
			iio_channel_attr_write_double(out0, "calibphase", (double) (-1.0 * vsin));
			iio_channel_attr_write_double(out1, "calibscale", (double) vcos);
			iio_channel_attr_write_double(out1, "calibphase", (double) vsin);
		}
	}
	return 0;
}

int streaming_interfaces(bool enable)
{
  if (enable)
  {
    rxa_chan_real = iio_device_find_channel(dev_rx, "voltage0", false);
    rxa_chan_imag = iio_device_find_channel(dev_rx, "voltage1", false);
    rxb_chan_real = iio_device_find_channel(dev_rx, "voltage4", false);
    rxb_chan_imag = iio_device_find_channel(dev_rx, "voltage5", false);
    if (!(rxa_chan_real && rxa_chan_imag && rxb_chan_real && rxb_chan_imag))
      streaming_interfaces(false);

    iio_channel_enable(rxa_chan_real);
    iio_channel_enable(rxa_chan_imag);
    iio_channel_enable(rxb_chan_real);
    iio_channel_enable(rxb_chan_imag);
    rxbuf = iio_device_create_buffer(dev_rx, SAMPLES, false);
    if (!rxbuf)
			streaming_interfaces(false);
  }
  else
  {
    if (rxbuf) { iio_buffer_destroy(rxbuf); }
    if (rxa_chan_real) { iio_channel_disable(rxa_chan_real); }
    if (rxa_chan_imag) { iio_channel_disable(rxa_chan_imag); }
    if (rxb_chan_real) { iio_channel_disable(rxb_chan_real); }
    if (rxb_chan_imag) { iio_channel_disable(rxb_chan_imag); }
		return -1;
  }
	return 0;
}

void read_buffer_data(struct iio_channel *chn, struct iio_buffer *buf, void *dst, size_t len)
{
  uintptr_t src_ptr, dst_ptr = (uintptr_t) dst, end = dst_ptr + len;
  unsigned int bytes = iio_channel_get_data_format(chn)->length / 8;
  uintptr_t buf_end = (uintptr_t) iio_buffer_end(buf);
  ptrdiff_t buf_step = iio_buffer_step(buf);

  for (src_ptr = (uintptr_t) iio_buffer_first(buf, chn);
    src_ptr < buf_end && dst_ptr + bytes <= end;
    src_ptr += buf_step, dst_ptr += bytes)
      iio_channel_convert(chn,(void *) dst_ptr, (const void *) src_ptr);

}

double estimate_phase_diff(double *estimate)
{
  ssize_t nbytes_rx = iio_buffer_refill(rxbuf);
  if (!nbytes_rx)
    return nbytes_rx;

  int16_t myData0_i[SAMPLES], myData0_q[SAMPLES];
  int16_t myData2_i[SAMPLES], myData2_q[SAMPLES];

  // Read data from all channels
  read_buffer_data(rxa_chan_real, rxbuf, myData0_i, SAMPLES*sizeof(int16_t));
  read_buffer_data(rxa_chan_imag, rxbuf, myData0_q, SAMPLES*sizeof(int16_t));
  read_buffer_data(rxb_chan_real, rxbuf, myData2_i, SAMPLES*sizeof(int16_t));
  read_buffer_data(rxb_chan_imag, rxbuf, myData2_q, SAMPLES*sizeof(int16_t));

  *estimate = calculate_phase(myData0_i,myData0_q,myData2_i,myData2_q, SAMPLES)*180/M_PI;
	return 0;
}

int calibrate_chain(struct iio_device *dev, double scale, double *phase)
{
	double est = 0, tmp;
  int k = 0, ret;
  for (;k<CALIBRATE_TRIES;k++)
  {
    *phase = STEP_SIZE*est+(*phase);
    ret = trx_phase_rotation(dev, *phase);
		CHECK(ret);

    if (streaming_interfaces(true)<0)
			return -ENODEV;
		ret = estimate_phase_diff(&est);
    CHECK(ret);
		est *= scale;
    streaming_interfaces(false);
    if (fabs(est)<TOLERANCE)
      break;
  }
#if (DEBUG > 0)
	printf("Remaining Phase error: %f\n",est);
  printf("Rotation: %f\n",*phase);
#endif

  return 0;
}

int quad_tracking(bool enable)
{
  struct iio_channel *chn = iio_device_find_channel(dev_phy, "voltage0", enable);
	if (chn == NULL)
		return -ENODEV;
  iio_channel_attr_write(chn, "quadrature_tracking_en", "0");
  chn = iio_device_find_channel(dev_phy_slave, "voltage0", enable);
	if (chn == NULL)
		return -ENODEV;
  iio_channel_attr_write(chn, "quadrature_tracking_en", "0");
	return 0;
}

int configure_transceiver(long long bw_hz,long long fs_hz,long long lo_hz)
{
  // Set up channels
  struct iio_channel *chnRX = iio_device_find_channel(dev_phy, "voltage0", false);
  struct iio_channel *chnTX = iio_device_find_channel(dev_phy, "voltage0", true);
	if (!(chnRX && chnRX))
		return -ENODEV;
  // Setup remaining channel info
	iio_channel_attr_write_longlong(chnRX, "rf_bandwidth", bw_hz);
	iio_channel_attr_write_longlong(chnRX, "sampling_frequency", fs_hz);
	iio_channel_attr_write_longlong(chnTX, "rf_bandwidth", bw_hz);
	iio_channel_attr_write_longlong(chnTX, "sampling_frequency", fs_hz);
  // Configure LO channel
  chnRX = iio_device_find_channel(dev_phy, "altvoltage0", true);
  chnTX = iio_device_find_channel(dev_phy, "altvoltage1", true);
	if (!(chnRX && chnRX))
		return -ENODEV;
	iio_channel_attr_write_longlong(chnRX, "frequency", lo_hz);
	iio_channel_attr_write_longlong(chnTX, "frequency", lo_hz);
	return 0;
}

int configure_dds(double fs, double scale)
{
  long long freq = (long long) fs*0.01;
	int i, j, ret = 0;

	for (i = 0; i < 2; i++) {
		for (j = 0; j < 8; j++) {
			ret |= iio_channel_attr_write_longlong(dds_out[i][j], "frequency", freq);
			ret |= iio_channel_attr_write_double(dds_out[i][j], "scale", scale);
		}

		dds_tx_phase_rotation(i ? dev_tx_slave : dev_tx, 0.0);
		trx_phase_rotation(i ? dev_tx_slave : dev_tx, 0.0);
	}
	return ret;
}

int get_dds_channels()
{
	struct iio_device *dev;
	int i, j;
	char name[16];

	for (i = 0; i < 2; i++) {
		dev = i ? dev_tx : dev_tx_slave;

		for (j = 0; j < 8; j++)
		{
			snprintf(name, sizeof(name), "altvoltage%d", j);

			dds_out[i][j] = iio_device_find_channel(dev, name, true);
			if (!dds_out[i][j])
			return -errno;
		}
	}
	return 0;
}

int setup_iio_devices(struct iio_context *ctx)
{
  dev_rx = iio_context_find_device(ctx, DEV_RX_NAME);
  dev_rx_slave = iio_context_find_device(ctx, DEV_RX_SLAVE_NAME);
  dev_phy = iio_context_find_device(ctx, DEV_PHY_NAME);
  dev_phy_slave = iio_context_find_device(ctx, DEV_PHY_SLAVE_NAME);
  dev_tx = iio_context_find_device(ctx, DEV_TX_NAME);
  dev_tx_slave = iio_context_find_device(ctx, DEV_TX_SLAVE_NAME);
	return (dev_rx&&dev_rx_slave&&dev_phy&&dev_phy_slave&&dev_tx&&dev_tx_slave);
}

/* Synchronize all transmit and receive channels for FMComms5*/
int phase_sync(struct iio_context *ctx, long long sample_rate, long long lo)
{
	// Set bandwidth same as sample rate
  long long bw = sample_rate;

  // Set up devices
  if (!setup_iio_devices(ctx))
		return -ENODEV;

  // Set up DDSs
  int ret = get_dds_channels();
	CHECK(ret);

  // Sync chips together
	ret = ad9361_multichip_sync(dev_phy, &dev_phy_slave, 1,
				FIXUP_INTERFACE_TIMING | CHECK_SAMPLE_RATES);
	CHECK(ret);

  // Set up DDS at given frequency
  ret = configure_dds(sample_rate, DDS_SCALE);
	CHECK(ret);

  // Set LO and bandwidths of transceivers
  ret = configure_transceiver(bw, sample_rate, lo);
	CHECK(ret);

	// Turn on quad tracking
  quad_tracking(true);

  // Align transmitter and receiver on chip A
  ret = trx_phase_rotation(dev_rx, 0.0);
	CHECK(ret);
	ret = trx_phase_rotation(dev_rx_slave, 0.0);
	CHECK(ret);
	configure_ports(1);
	double phase_est_last=0, phase_est=0;
  ret = calibrate_chain(dev_rx_slave, -1, &phase_est_last);
	CHECK(ret);

  // Align receivers across chip
  ret = trx_phase_rotation(dev_rx_slave, 0.0);
	CHECK(ret);
  configure_ports(3);
  ret = calibrate_chain(dev_rx, 1, &phase_est);
	CHECK(ret);
  // Align transmitters
  ret = trx_phase_rotation(dev_rx_slave, 0);
	CHECK(ret);
  configure_ports(4);
  ret = calibrate_chain(dev_tx_slave, -1, &phase_est);
	CHECK(ret);

  // Set rotation of chip B receiver
  ret = trx_phase_rotation(dev_rx_slave, phase_est_last);
	CHECK(ret);

  return 0;
}

/* Synchronize all transmit and receive channels for FMComms5*/
int ad9361_fmcomms5_phase_sync(struct iio_context *ctx, long long sample_rate, long long lo)
{
	int ret = phase_sync(ctx, sample_rate, lo);
	// Reset ports out to RF
  configure_ports(0);
	return ret;
}
