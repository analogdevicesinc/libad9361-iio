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

#include "ad9361.h"

#include <errno.h>
#include <iio.h>
#include <stdio.h>
#ifdef _WIN32
#include <windows.h>
#else
#include <unistd.h>
#include <time.h>
#endif

#define MAX_AD9361_SYNC_DEVS	4

static void ad9361_sleep_ms(void)
{
#ifdef _WIN32
	Sleep(1); /* milliseconds */
#else
	struct timespec time;

	time.tv_sec = 0;
	time.tv_nsec = 1000 * 1000;
	nanosleep(&time, NULL);
#endif
}

int ad9361_multichip_sync(struct iio_device *master, struct iio_device **slaves,
		unsigned int num_slaves, unsigned int flags)
{
	char ensm_mode[MAX_AD9361_SYNC_DEVS][20];
	unsigned int i, step;
	bool mcs_is_debug_attr = !iio_device_find_attr(master, "multichip_sync");

	if (num_slaves >= MAX_AD9361_SYNC_DEVS || num_slaves < 1)
		return -EINVAL;

	if (flags & CHECK_SAMPLE_RATES) {
		struct iio_channel *tx_sample_master, *tx_sample_slave;
		long long tx_sample_master_freq, tx_sample_slave_freq;

		tx_sample_master = iio_device_find_channel(master, "voltage0", true);
		iio_channel_attr_read_longlong(tx_sample_master, "sampling_frequency", &tx_sample_master_freq);

		for (i = 0; i < num_slaves; i++) {
			tx_sample_slave = iio_device_find_channel(slaves[i], "voltage0", true);
			if (tx_sample_slave == NULL)
				return -ENODEV;

			iio_channel_attr_read_longlong(tx_sample_slave, "sampling_frequency", &tx_sample_slave_freq);

			if (tx_sample_master_freq != tx_sample_slave_freq) {
				fprintf(stderr, "tx_sample_master_freq != tx_sample_slave_freq\nUpdating...\n");
				iio_channel_attr_write_longlong(tx_sample_slave, "sampling_frequency", tx_sample_master_freq);
			}
		}
	}

	if (flags & FIXUP_INTERFACE_TIMING) {
		unsigned tmp, tmp2;
		iio_device_reg_read(master, 0x6, &tmp);
		iio_device_reg_read(master, 0x7, &tmp2);

		for (i = 0; i < num_slaves; i++) {
			iio_device_reg_write(slaves[i], 0x6, tmp);
			iio_device_reg_write(slaves[i], 0x7, tmp2);
		}
	}

	/* Move the parts int ALERT for MCS */
	iio_device_attr_read(master, "ensm_mode", ensm_mode[0], sizeof(ensm_mode));
	iio_device_attr_write(master, "ensm_mode", "alert");

	for (i = 0; i < num_slaves; i++) {
		iio_device_attr_read(slaves[i], "ensm_mode", ensm_mode[i + 1], sizeof(ensm_mode));
		iio_device_attr_write(slaves[i], "ensm_mode", "alert");
	}

	for (step = 0; step <= 5; step++) {
		for (i = 0; i < num_slaves; i++) {
			if (mcs_is_debug_attr)
				iio_device_debug_attr_write_longlong(slaves[i], "multichip_sync", step);
			else
				iio_device_attr_write_longlong(slaves[i], "multichip_sync", step);
		}

		/* The master controls the SYNC GPIO */
		if (mcs_is_debug_attr)
			iio_device_debug_attr_write_longlong(master, "multichip_sync", step);
		else
			iio_device_attr_write_longlong(master, "multichip_sync", step);

		ad9361_sleep_ms();
	}

	iio_device_attr_write(master, "ensm_mode", ensm_mode[0]);

	for (i = 0; i < num_slaves; i++)
		iio_device_attr_write(slaves[i], "ensm_mode", ensm_mode[i + 1]);

	return 0;
}

int ad9361_fmcomms5_multichip_sync(struct iio_context *ctx, unsigned int flags)
{
	struct iio_device *master, *slave;

	master = iio_context_find_device(ctx, "ad9361-phy");
	slave = iio_context_find_device(ctx, "ad9361-phy-B");

	if (!master || !slave)
		return -ENODEV;

	return ad9361_multichip_sync(master, &slave, 1, flags);
}
