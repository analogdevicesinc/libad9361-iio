#!/usr/bin/env python
# SPDX-License-Identifier: LGPL-2.1-or-later
"""
SPDX-License-Identifier: LGPL-2.1-or-later
Copyright (C) 2021 Analog Devices, Inc.
Author: Travis F. Collins <travis.collins@analog.com>
"""

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


def multichip_sync(main: iio.Device, secondaries: List[iio.Device], flags: int):
    """ """

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
    """ """
    if math.trunc(flags) != flags:
        raise Exception("flags must be an integer")

    cflags = c_uint(flags)

    ret = _ad9361_fmcomms5_multichip_sync(ctx._context, cflags)

    if ret != 0:
        raise Exception(f"Setting MCS failed. ERROR: {ret}")


def set_bb_rate(dev: iio.Device, rate: int):
    """ """
    if math.trunc(rate) != rate:
        raise Exception("rate must be an integer")

    ret = _ad9361_set_bb_rate(dev._device, c_ulong(int(rate)))

    if ret != 0:
        raise Exception(f"Setting rate failed. ERROR: {ret}")
