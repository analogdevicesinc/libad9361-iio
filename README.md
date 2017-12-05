# libad9361-iio

[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/analogdevicesinc/libad9361-iio?branch=master&svg=true)](https://ci.appveyor.com/project/analogdevicesinc/libad9361-iio)
[![Travis Build Status](https://travis-ci.org/analogdevicesinc/libad9361-iio.svg?branch=master)](https://travis-ci.org/analogdevicesinc/libad9361-iio)
[![License](https://img.shields.io/badge/License-LGPL%20v2.1-blue.svg)](https://github.com/analogdevicesinc/libad9361-iio/blob/master/LICENSE)

This is a simple library used for userspace, which manages multi-chip sync,
on platforms (FMCOMMS5) where multiple AD9361 devices are used.

## Building & Installing

should be a quick matter of `cmake`, then `make`, then `make install`:

```
rgetz@pinky:~/libad9361-iio$ cmake ./CMakeLists.txt
-- The C compiler identification is GNU 4.7.2
-- Check for working C compiler: /usr/bin/gcc
-- Check for working C compiler: /usr/bin/gcc -- works
-- Detecting C compiler ABI info
-- Detecting C compiler ABI info - done
-- Configuring done
-- Generating done
-- Build files have been written to: /home/rgetz/libad9361-iio
rgetz@pinky:~/libad9361-iio$ make
Scanning dependencies of target ad9361
[100%] Building C object CMakeFiles/ad9361.dir/ad9361_multichip_sync.c.o
Linking C shared library libad9361.so
Copying OS X content Headers/ad9361.h
[100%] Built target ad9361
rgetz@pinky:~/libad9361-iio$ sudo make install
[sudo] password for rgetz: 
[100%] Built target ad9361
Install the project...
-- Install configuration: ""
-- Installing: /usr/local/lib/pkgconfig/libad9361.pc
-- Installing: /usr/local/lib/libad9361.so.0.1
-- Installing: /usr/local/lib/libad9361.so.0
-- Installing: /usr/local/lib/libad9361.so
-- Installing: /usr/local/include/ad9361.h
```
