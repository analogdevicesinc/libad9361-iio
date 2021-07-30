#!/usr/bin/env python
# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright (C) 2021 Analog Devices, Inc.
# Author: Travis F. Collins <travis.collins@analog.com>

from typing import List
import math

import iio
from iio import _ContextPtr, _DevicePtr

# Imports from package ctypes are not grouped
# pylint: disable=ungrouped-imports
#
# The way the methods and classes are used, violate this, but there
#   isn't a good way to do things otherwise
# pylint: disable=protected-access
from ctypes import (
    c_longlong,
    c_uint,
    c_int,
    c_ulong,
    POINTER as _POINTER,
    CDLL as _cdll,
)
from ctypes.util import find_library
from os import strerror as _strerror
from platform import system as _system

if "Windows" in _system():
    from ctypes import get_last_error
else:
    from ctypes import get_errno
# pylint: enable=ungrouped-imports

# ctypes requires the errcheck to take three arguments, but we don't use them
# pylint: disable=unused-argument


def _check_negative(result, func, arguments):
    if result >= 0:
        return result
    raise OSError(-result, _strerror(-result))


if "Windows" in _system():
    _libad9361 = "libad9361.dll"
else:
    # Non-windows, possibly Posix system
    _libad9361 = "ad9361"

_lib = _cdll(find_library(_libad9361), use_errno=True, use_last_error=True)


_ad9361_multichip_sync = _lib.ad9361_multichip_sync
_ad9361_multichip_sync.restype = c_int
_ad9361_multichip_sync.argtypes = (_DevicePtr, _POINTER(_DevicePtr), c_uint, c_uint)
_ad9361_multichip_sync.errcheck = _check_negative

_ad9361_fmcomms5_multichip_sync = _lib.ad9361_fmcomms5_multichip_sync
_ad9361_fmcomms5_multichip_sync.restype = c_int
_ad9361_fmcomms5_multichip_sync.argtypes = (_ContextPtr, c_uint)
_ad9361_fmcomms5_multichip_sync.errcheck = _check_negative

_ad9361_set_bb_rate = _lib.ad9361_set_bb_rate
_ad9361_set_bb_rate.restype = c_int
_ad9361_set_bb_rate.argtypes = (_DevicePtr, c_ulong)
_ad9361_set_bb_rate.errcheck = _check_negative

_ad9361_set_trx_fir_enable = _lib.ad9361_set_trx_fir_enable
_ad9361_set_trx_fir_enable.restype = c_int
_ad9361_set_trx_fir_enable.argtypes = (_DevicePtr, c_int)
_ad9361_set_trx_fir_enable.errcheck = _check_negative

_ad9361_get_trx_fir_enable = _lib.ad9361_get_trx_fir_enable
_ad9361_get_trx_fir_enable.restype = c_int
_ad9361_get_trx_fir_enable.argtypes = (_DevicePtr, c_int)
_ad9361_get_trx_fir_enable.errcheck = _check_negative

# _ad9361_generate_fir_taps = _lib.ad9361_generate_fir_taps
# _ad9361_generate_fir_taps.restype = c_int
# _ad9361_generate_fir_taps.argtypes = (_DevicePtr, c_int)
# _ad9361_generate_fir_taps.errcheck = _check_negative

# _ad9361_calculate_rf_clock_chain = _lib.ad9361_calculate_rf_clock_chain
# _ad9361_calculate_rf_clock_chain.restype = c_int
# _ad9361_calculate_rf_clock_chain.argtypes = (_DevicePtr, c_int)
# _ad9361_calculate_rf_clock_chain.errcheck = _check_negative

# _ad9361_calculate_rf_clock_chain_fdp = _lib.ad9361_calculate_rf_clock_chain_fdp
# _ad9361_calculate_rf_clock_chain_fdp.restype = c_int
# _ad9361_calculate_rf_clock_chain_fdp.argtypes = (_DevicePtr, c_int)
# _ad9361_calculate_rf_clock_chain_fdp.errcheck = _check_negative

_ad9361_set_bb_rate_custom_filter_auto = _lib.ad9361_set_bb_rate_custom_filter_auto
_ad9361_set_bb_rate_custom_filter_auto.restype = c_int
_ad9361_set_bb_rate_custom_filter_auto.argtypes = (_DevicePtr, c_ulong)
_ad9361_set_bb_rate_custom_filter_auto.errcheck = _check_negative

# _ad9361_set_bb_rate_custom_filter_manual = _lib.ad9361_set_bb_rate_custom_filter_manual
# _ad9361_set_bb_rate_custom_filter_manual.restype = c_int
# _ad9361_set_bb_rate_custom_filter_manual.argtypes = (_DevicePtr, c_int)
# _ad9361_set_bb_rate_custom_filter_manual.errcheck = _check_negative

_ad9361_fmcomms5_phase_sync = _lib.ad9361_fmcomms5_phase_sync
_ad9361_fmcomms5_phase_sync.restype = c_int
_ad9361_fmcomms5_phase_sync.argtypes = (_ContextPtr, c_longlong)
_ad9361_fmcomms5_phase_sync.errcheck = _check_negative


def multichip_sync(main: iio.Device, secondaries: List[iio.Device], flags: int) -> None:
    """Multi-chip synchronization (MCS) management.

    :param: iio.Device master: IIO Device object of master AD9361 transceiver
    :param: List[iio.Device] slaves: A list of IIO Device objects to child AD9361 transceivers
    :param int flags: Control flags for MCS configuration
    """

    if math.trunc(flags) != flags:
        raise Exception("flags must be an integer")

    psecondaries = [d._device for d in secondaries]
    csecondaries = (_DevicePtr * len(secondaries))(*psecondaries)

    cnum_secondaries = c_uint(len(secondaries))
    cflags = c_uint(flags)

    ret = _ad9361_multichip_sync(main._device, csecondaries, cnum_secondaries, cflags)

    if ret != 0:
        raise Exception(f"Setting MCS failed. ERROR: {ret}")


def fmcomms5_multichip_sync(ctx: iio.Context, flags: int):
    """FMComms5 Multi-chip synchronization (MCS).

    :param iio.Context ctx: IIO Context object of FMComms5
    :param int flags: Control flags for MCS configuration
    """
    if math.trunc(flags) != flags:
        raise Exception("flags must be an integer")

    cflags = c_uint(flags)

    ret = _ad9361_fmcomms5_multichip_sync(ctx._context, cflags)

    if ret != 0:
        raise Exception(f"Setting MCS failed. ERROR: {ret}")


def set_bb_rate(dev: iio.Device, rate: int):
    """Baseband rate configuration with generic filter support.

    :param: iio.Device dev: IIO Device object of AD9361 transceiver phy driver
    :param int rate: Desired sample rate between 61.44 MSPS and 520833 SPS
    """
    if math.trunc(rate) != rate:
        raise Exception("rate must be an integer")

    ret = _ad9361_set_bb_rate(dev._device, c_ulong(int(rate)))

    if ret != 0:
        raise Exception(f"Setting rate failed. ERROR: {ret}")


def set_trx_fir_enable(dev: iio.Device, enable: int):
    """Enable or disable transmit and receiver FIRs simultaneously.

    :param: iio.Device dev: IIO Device object of AD9361 transceiver phy driver
    :param int enable: Enable FIRs when 1 or disable when 0
    """
    if math.trunc(enable) != enable:
        raise Exception("enable must be an integer")

    ret = _ad9361_set_trx_fir_enable(dev._device, c_int(int(enable)))

    if ret != 0:
        raise Exception(f"Failed setting filters state. ERROR: {ret}")


def get_trx_fir_enable(dev: iio.Device):
    """Get current enable value of transmit and receiver FIRs.

    :param: iio.Device dev: IIO Device object of AD9361 transceiver phy driver
    """
    enable = _POINTER(c_int(0))
    ret = _ad9361_set_trx_fir_enable(dev._device, enable)

    if ret != 0:
        raise Exception(f"Failed to get filters state. ERROR: {ret}")

    return int(enable)


def set_bb_rate_custom_filter_auto(dev: iio.Device, rate: int):
    """Baseband rate configuration with custom filter support based on desired
    baseband sample rate.

    :param: iio.Device dev: IIO Device object of AD9361 transceiver phy driver
    :param int rate: Desired sample rate between 61.44 MSPS and 520833 SPS
    """
    if math.trunc(rate) != rate:
        raise Exception("rate must be an integer")

    ret = _ad9361_set_bb_rate_custom_filter_auto(dev._device, c_ulong(int(rate)))

    if ret != 0:
        raise Exception(f"Failed to set custom filter. ERROR: {ret}")


def fmcomms5_phase_sync(ctx: iio.Context, lo: int):
    """FMComms5 phase synchronize all TX and RX channels together.

    :param iio.Context ctx: IIO Context object of FMComms5
    :param int lo: Frequency in hertz of LO for TX and RX
    """
    if math.trunc(lo) != lo:
        raise Exception("lo must be an integer")

    ret = _ad9361_fmcomms5_phase_sync(ctx._ctx, c_longlong(int(lo)))

    if ret != 0:
        raise Exception(f"Failed to sync FMComms5. ERROR: {ret}")
