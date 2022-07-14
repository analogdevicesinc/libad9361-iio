# libad9361-iio: Python Bindings

This package contains the python bindings for libad9361-iio, a library for managing different AD9361 transceiver features. Such as multi-chip sync and filter generation.

License : [![License](https://img.shields.io/badge/license-LGPL2-blue.svg)](https://github.com/analogdevicesinc/libad9361-iio/blob/master/COPYING.txt)
Latest Release : [![GitHub release](https://img.shields.io/github/release/analogdevicesinc/libad9361-iio.svg)](https://github.com/analogdevicesinc/libad9361-iio/releases/latest)
Downloads :  [![Github All Releases](https://img.shields.io/github/downloads/analogdevicesinc/libad9361-iio/total.svg)](https://github.com/analogdevicesinc/libad9361-iio/releases/latest)

Support:<br>
If you have a question about libad9361-iio or the AD9361 transceiver please ask on : [![EngineerZone](https://img.shields.io/badge/chat-on%20EngineerZone-blue.svg)](https://ez.analog.com/linux-device-drivers/linux-software-drivers).

## Requirements
To use these bindings naturally you need the core library they depend upon, libad9361-iio. This is not packaged with the pypi release but there are a number of options:
  - If you want to just use libad9361-iio, we suggest using the [latest release](https://github.com/analogdevicesinc/libad9361-iio/releases/latest).
  - If you think you have found a bug in the release, or need a feature which isn't in the release, try the latest **untested** binaries from the master branch and check out the [documentation](https://codedocs.xyz/analogdevicesinc/libad9361-iio/) based on the master branch. We provide builds for a few operating systems. If you need something else, we can most likely add that -- just ask.

### Installing the bindings
To install these bindings there are a few methods. If you already have the library itself and just need the bindings, pip is the most convenient method:
```shell
(sudo) pip install pylibad9361-iio
```
If you do not want to use pip, then installation is dependent on your operating system.
#### Linux / macOS
For Linux and macOS the python bindings need to be installed through source if not using pip.

#### Windows
Only pip installation is supported.

### Support
If you have a question about libad9361-iio or the python bindings please ask on: [![EngineerZone](https://img.shields.io/badge/chat-on%20EngineerZone-blue.svg)](https://ez.analog.com/linux-device-drivers/linux-software-drivers).

If you use it, and like it - please let us know. If you use it, and hate it - please let us know that too. The goal of the project is to try to make the AD9361 transceiver easier to use on a variety of platforms. If we aren't doing that - we will try to make it better.

Feedback is appreciated (in order of preference):

  * [Github trackers](https://github.com/analogdevicesinc/libad9361-iio/issues) for bugs, improvements, or feature requests
  * [Analog Devices web forums](https://ez.analog.com/community/linux-device-drivers/linux-software-drivers) for general help on libiio and/or ADI Linux IIO drivers

## Useful resources
  * [API Documentation](http://analogdevicesinc.github.io/libad9361-iio/)
  * [Libiio](http://wiki.analog.com/resources/tools-software/linux-software/libiio)

