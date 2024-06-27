## :warning: Important note

Currently the main branch of this repository does not support the new v1.0 API of [libiio](https://github.com/analogdevicesinc/libiio).

If you want to build libad9361-iio please use the [libad9361-iio-v0](https://github.com/analogdevicesinc/libad9361-iio/tree/libad9361-iio-v0) branch that supports the old v0.x API, that can be found on the [libiio-v0](https://github.com/analogdevicesinc/libiio/tree/libiio-v0) branch or as packages on the release section of the libiio repository.


# libad9361-iio

This is a simple library used for userspace,
 - which manages multi-chip sync (on platforms (FMCOMMS5) where multiple AD9361 devices are use)
 - can create AD9361 specific FIR filters on the fly,

**Dependencies**

 - [libiio IIO interface library (Version less than v1.0)](https://github.com/analogdevicesinc/libiio)

**Docs**

Doxygen-based documentation is available at: http://analogdevicesinc.github.io/libad9361-iio/


License : [![License](https://img.shields.io/badge/license-LGPL2-blue.svg)](https://github.com/analogdevicesinc/libad9361-iio/blob/master/COPYING.txt)
Latest Release : [![GitHub release](https://img.shields.io/github/release/analogdevicesinc/libad9361-iio.svg)](https://github.com/analogdevicesinc/libad9361-iio/releases/latest)
Downloads :  [![Github All Releases](https://img.shields.io/github/downloads/analogdevicesinc/libad9361-iio/total.svg)](https://github.com/analogdevicesinc/libad9361-iio/releases/latest)

As with many open source packages, we use [GitHub](https://github.com/analogdevicesinc/libad9361-iio) to do develop and maintain the source, and [Azure Pipelines](https://azure.microsoft.com/en-gb/services/devops/pipelines/) for continuous integration.
  - If you want to just use libad9361-iio, we suggest using the [latest release](https://github.com/analogdevicesinc/libad9361-iio/releases/latest).
  - If you think you have found a bug in the release, or need a feature which isn't in the release, try the latest **untested** binaries from the master branch. We provide builds for a few operating systems. If you need something else, we can most likely add that -- just ask.


If you use it, and like it - please let us know. If you use it, and hate it - please let us know that too. The goal of the project is to try to make Linux IIO devices easier to use on a variety of platforms. If we aren't doing that - we will try to make it better.


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
