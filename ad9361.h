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

/* FLAGS */
#define FIXUP_INTERFACE_TIMING	1
#define CHECK_SAMPLE_RATES	2

struct iio_context;
struct iio_device;

int ad9361_multichip_sync(struct iio_device *master, struct iio_device **slaves,
		unsigned int num_slaves, unsigned int flags);

int ad9361_fmcomms5_multichip_sync(struct iio_context *ctx, unsigned int flags);

#endif /* __AD9361_H__ */
