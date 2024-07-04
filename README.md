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

| Operating System        | GitHub libad9361-ii-v0 status  | Version |  Installer Package  | tarball or zip |
|:-----------------------:|:---------------------:|:-------:|:-------------------:|:--------------:|
| Windows                 | [![Build Status](https://dev.azure.com/AnalogDevices/OpenSource/_apis/build/status%2Fanalogdevicesinc.libad9361-iio?branchName=libad9361-iio-v0&jobName=WindowsBuilds&configuration=WindowsBuilds%20VS2019_Win64)](https://dev.azure.com/AnalogDevices/OpenSource/_build/latest?definitionId=10&branchName=libad9361-iio-v0) | Windows-64 Server 2019 | [![Latest Windows installer](https://raw.githubusercontent.com/wiki/analogdevicesinc/libiio/img/win_box.png)](https://swdownloads.analog.com/cse/azure_builds/libad9361-setup.exe) | [![Latest 64-bit Windows zip](https://raw.githubusercontent.com/wiki/analogdevicesinc/libiio/img/win_box.png)](https://swdownloads.analog.com/cse/azure_builds/Windows-VS-2019-x64-latest_libad9361-iio-v0_libad9361.zip) |
|                         | [![Build Status](https://dev.azure.com/AnalogDevices/OpenSource/_apis/build/status%2Fanalogdevicesinc.libad9361-iio?branchName=libad9361-iio-v0&jobName=WindowsBuilds&configuration=WindowsBuilds%20VS2022)](https://dev.azure.com/AnalogDevices/OpenSource/_build/latest?definitionId=10&branchName=libad9361-iio-v0) | Windows-64 Server 2022 | [![Latest Windows installer](https://raw.githubusercontent.com/wiki/analogdevicesinc/libiio/img/win_box.png)](https://swdownloads.analog.com/cse/azure_builds/libad9361-setup.exe) | [![Latest 64-bit Windows zip](https://raw.githubusercontent.com/wiki/analogdevicesinc/libiio/img/win_box.png)](https://swdownloads.analog.com/cse/azure_builds/Windows-VS-2022-x64-latest_libad9361-iio-v0_libad9361.zip) |
| OS X                    | [![Build Status](https://dev.azure.com/AnalogDevices/OpenSource/_apis/build/status%2Fanalogdevicesinc.libad9361-iio?branchName=libad9361-iio-v0&jobName=macOSBuilds&configuration=macOSBuilds%20macOS_13_x64)](https://dev.azure.com/AnalogDevices/OpenSource/_build/latest?definitionId=10&branchName=libad9361-iio-v0) | macOS Ventura <br />(v 13 x64) | [![OS-X package 10.14](https://raw.githubusercontent.com/wiki/analogdevicesinc/libiio/img/osx_box.png)](https://swdownloads.analog.com/cse/azure_builds/macOS-13-x64_latest_libad9361-iio-v0_libad9361.pkg) | |
|                         | [![Build Status](https://dev.azure.com/AnalogDevices/OpenSource/_apis/build/status%2Fanalogdevicesinc.libad9361-iio?branchName=libad9361-iio-v0&jobName=macOSBuilds&configuration=macOSBuilds%20macOS_13_arm64)](https://dev.azure.com/AnalogDevices/OpenSource/_build/latest?definitionId=10&branchName=libad9361-iio-v0) |  macOS Ventura <br />(v 13 arm64) | [![OS-X package 10.13](https://raw.githubusercontent.com/wiki/analogdevicesinc/libiio/img/osx_box.png)](https://swdownloads.analog.com/cse/azure_builds/macOS-13-arm64_latest_libad9361-iio-v0_libad9361.pkg) | |
|                    | [![Build Status](https://dev.azure.com/AnalogDevices/OpenSource/_apis/build/status%2Fanalogdevicesinc.libad9361-iio?branchName=libad9361-iio-v0&jobName=macOSBuilds&configuration=macOSBuilds%20macOS_12)](https://dev.azure.com/AnalogDevices/OpenSource/_build/latest?definitionId=10&branchName=libad9361-iio-v0) | macOS Monterey <br />(v 12) | [![OS-X package 10.12](https://raw.githubusercontent.com/wiki/analogdevicesinc/libiio/img/osx_box.png)](https://swdownloads.analog.com/cse/azure_builds/macOS-12_latest_libad9361-iio-v0_libad9361.pkg) | |
| Linux     | [![Build Status](https://dev.azure.com/AnalogDevices/OpenSource/_apis/build/status%2Fanalogdevicesinc.libad9361-iio?branchName=libad9361-iio-v0&jobName=LinuxBuilds&configuration=LinuxBuilds%20ubuntu_22_04_x86_64)](https://dev.azure.com/AnalogDevices/OpenSource/_build/latest?definitionId=10&branchName=libad9361-iio-v0) |  Ubuntu Jammy Jellyfish<br />(v 22.04)<sup>1</sup>  | [![Debian](https://raw.githubusercontent.com/wiki/analogdevicesinc/libiio/img/deb.png)](https://swdownloads.analog.com/cse/azure_builds/Ubuntu-22.04_latest_libad9361-iio-v0_libad9361.deb) | |
|  | [![Build Status](https://dev.azure.com/AnalogDevices/OpenSource/_apis/build/status%2Fanalogdevicesinc.libad9361-iio?branchName=libad9361-iio-v0&jobName=LinuxBuilds&configuration=LinuxBuilds%20ubuntu_20_04_x86_64)](https://dev.azure.com/AnalogDevices/OpenSource/_build/latest?definitionId=10&branchName=libad9361-iio-v0) | Ubuntu Focal Fossa<br />(v 20.04)<sup>1</sup> | [![Debian](https://raw.githubusercontent.com/wiki/analogdevicesinc/libiio/img/deb.png)](https://swdownloads.analog.com/cse/azure_builds/Ubuntu-20.04_latest_libad9361-iio-v0_libad9361.deb) | |
|  | [![Build Status](https://dev.azure.com/AnalogDevices/OpenSource/_apis/build/status%2Fanalogdevicesinc.libad9361-iio?branchName=libad9361-iio-v0&jobName=LinuxBuilds&configuration=LinuxBuilds%20fedora34)](https://dev.azure.com/AnalogDevices/OpenSource/_build/latest?definitionId=10&branchName=libad9361-iio-v0) | Fedora 34  | [![RPM File](https://raw.githubusercontent.com/wiki/analogdevicesinc/libiio/img/rpm.png)](https://swdownloads.analog.com/cse/azure_builds/Fedora-34_latest_libad9361-iio-v0_libad9361.rpm)  | |
|  | [![Build Status](https://dev.azure.com/AnalogDevices/OpenSource/_apis/build/status%2Fanalogdevicesinc.libad9361-iio?branchName=libad9361-iio-v0&jobName=LinuxBuilds&configuration=LinuxBuilds%20fedora28)](https://dev.azure.com/AnalogDevices/OpenSource/_build/latest?definitionId=10&branchName=libad9361-iio-v0) | Fedora 28  | [![RPM File](https://raw.githubusercontent.com/wiki/analogdevicesinc/libiio/img/rpm.png)](https://swdownloads.analog.com/cse/azure_builds/Fedora-28_latest_libad9361-iio-v0_libad9361.rpm)  | |
|  | [![Build Status](https://dev.azure.com/AnalogDevices/OpenSource/_apis/build/status%2Fanalogdevicesinc.libad9361-iio?branchName=libad9361-iio-v0&jobName=LinuxBuilds&configuration=LinuxBuilds%20centos_7)](https://dev.azure.com/AnalogDevices/OpenSource/_build/latest?definitionId=10&branchName=libad9361-iio-v0) | CentOS 7  | [![RPM File](https://raw.githubusercontent.com/wiki/analogdevicesinc/libiio/img/rpm.png)](https://swdownloads.analog.com/cse/azure_builds/CentOS-7_latest_libad9361-iio-v0_libad9361.rpm)  |  |
|  | [![Build Status](https://dev.azure.com/AnalogDevices/OpenSource/_apis/build/status%2Fanalogdevicesinc.libad9361-iio?branchName=libad9361-iio-v0&jobName=LinuxBuilds&configuration=LinuxBuilds%20debian_bullseye)](https://dev.azure.com/AnalogDevices/OpenSource/_build/latest?definitionId=10&branchName=libad9361-iio-v0) | Debian Bullseye | [![Debian](https://raw.githubusercontent.com/wiki/analogdevicesinc/libiio/img/deb.png)](https://swdownloads.analog.com/cse/azure_builds/Debian-11_latest_libad9361-iio-v0_libad9361.deb) |  |
|  | [![Build Status](https://dev.azure.com/AnalogDevices/OpenSource/_apis/build/status%2Fanalogdevicesinc.libad9361-iio?branchName=libad9361-iio-v0&jobName=LinuxBuilds&configuration=LinuxBuilds%20debian_bookworm)](https://dev.azure.com/AnalogDevices/OpenSource/_build/latest?definitionId=10&branchName=libad9361-iio-v0) | Debian Bookworm | [![Debian](https://raw.githubusercontent.com/wiki/analogdevicesinc/libiio/img/deb.png)](https://swdownloads.analog.com/cse/azure_builds/Debian-12_latest_libad9361-iio-v0_libad9361.deb) |  |
|  | [![Build Status](https://dev.azure.com/AnalogDevices/OpenSource/_apis/build/status%2Fanalogdevicesinc.libad9361-iio?branchName=libad9361-iio-v0&jobName=LinuxBuilds&configuration=LinuxBuilds%20opensuse_15_4)](https://dev.azure.com/AnalogDevices/OpenSource/_build/latest?definitionId=10&branchName=libad9361-iio-v0) | openSUSE 15.4 | [![Debian](https://raw.githubusercontent.com/wiki/analogdevicesinc/libiio/img/deb.png)](https://swdownloads.analog.com/cse/azure_builds/openSUSE-15.4_latest_libad9361-iio-v0_libad9361.deb) |  |
| ARM     | [![Build Status](https://dev.azure.com/AnalogDevices/OpenSource/_apis/build/status%2Fanalogdevicesinc.libad9361-iio?branchName=libad9361-iio-v0&jobName=ARMBuilds&configuration=ARMBuilds%20ubuntu-ppc64le)](https://dev.azure.com/AnalogDevices/OpenSource/_build/latest?definitionId=10&branchName=libad9361-iio-v0) | Ubuntu-ppc64le | [![Debian](https://raw.githubusercontent.com/wiki/analogdevicesinc/libiio/img/deb.png)](https://swdownloads.analog.com/cse/azure_builds/Ubuntu-ppc64le_latest_libad9361-iio-v0_libad9361.deb) |  |
|  | [![Build Status](https://dev.azure.com/AnalogDevices/OpenSource/_apis/build/status%2Fanalogdevicesinc.libad9361-iio?branchName=libad9361-iio-v0&jobName=ARMBuilds&configuration=ARMBuilds%20ubuntu-x390x)](https://dev.azure.com/AnalogDevices/OpenSource/_build/latest?definitionId=10&branchName=libad9361-iio-v0) | Ubuntu-x390x | [![Debian](https://raw.githubusercontent.com/wiki/analogdevicesinc/libiio/img/deb.png)](https://swdownloads.analog.com/cse/azure_builds/Ubuntu-x390x_latest_libad9361-iio-v0_libad9361.deb) |  |
|  | [![Build Status](https://dev.azure.com/AnalogDevices/OpenSource/_apis/build/status%2Fanalogdevicesinc.libad9361-iio?branchName=libad9361-iio-v0&jobName=ARMBuilds&configuration=ARMBuilds%20debian_buster_arm64v8)](https://dev.azure.com/AnalogDevices/OpenSource/_build/latest?definitionId=10&branchName=libad9361-iio-v0) | Ubuntu-arm64v8 | [![Debian](https://raw.githubusercontent.com/wiki/analogdevicesinc/libiio/img/deb.png)](https://swdownloads.analog.com/cse/azure_builds/Ubuntu-arm64v8_latest_libad9361-iio-v0_libad9361.deb) |  |
|  | [![Build Status](https://dev.azure.com/AnalogDevices/OpenSource/_apis/build/status%2Fanalogdevicesinc.libad9361-iio?branchName=libad9361-iio-v0&jobName=ARMBuilds&configuration=ARMBuilds%20debian_buster_arm32v7)](https://dev.azure.com/AnalogDevices/OpenSource/_build/latest?definitionId=10&branchName=libad9361-iio-v0) | Ubuntu-arm32v7 | [![Debian](https://raw.githubusercontent.com/wiki/analogdevicesinc/libiio/img/deb.png)](https://swdownloads.analog.com/cse/azure_builds/Ubuntu-arm32v7_latest_libad9361-iio-v0_libad9361.deb) |  |


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
