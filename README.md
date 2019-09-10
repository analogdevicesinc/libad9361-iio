# libad9361-iio

This is a simple library used for userspace,
 - which manages multi-chip sync (on platforms (FMCOMMS5) where multiple AD9361 devices are use)
 - can create AD9361 specific FIR filters on the fly,
 
**Docs**

Doxygen-based documentation is available at: http://analogdevicesinc.github.io/libad9361-iio/


License : [![License](https://img.shields.io/badge/license-LGPL2-blue.svg)](https://github.com/analogdevicesinc/libad9361-iio/blob/master/COPYING.txt)
Latest Release : [![GitHub release](https://img.shields.io/github/release/analogdevicesinc/libad9361-iio.svg)](https://github.com/analogdevicesinc/libad9361-iio/releases/latest)
Downloads :  [![Github All Releases](https://img.shields.io/github/downloads/analogdevicesinc/libad9361-iio/total.svg)](https://github.com/analogdevicesinc/libad9361-iio/releases/latest)

As with many open source packages, we use [GitHub](https://github.com/analogdevicesinc/libad9361-iio) to do develop and maintain the source, and [Travis CI](https://travis-ci.com/) and [Appveyor](https://www.appveyor.com/) for continuous integration.
  - If you want to just use libad9361-iio, we suggest using the [latest release](https://github.com/analogdevicesinc/libad9361-iio/releases/latest).
  - If you think you have found a bug in the release, or need a feature which isn't in the release, try the latest **untested** binaries from the master branch. We provide builds for a few operating systems. If you need something else, we can most likely add that -- just ask.

| Operating System        | GitHub master status  | Version |  Installer Package  | tarball or zip |
|:-----------------------:|:---------------------:|:-------:|:-------------------:|:--------------:|
| Windows                 | [![Windows Status](https://ci.appveyor.com/api/projects/status/github/analogdevicesinc/libad9361-iio?svg=true)](https://ci.appveyor.com/project/analogdevicesinc/libad9361-iio/branch/master) | Windows 10<br />Windows 8.1<br />Windows 8<br />Windows 7 | [![Latest Windows installer](https://raw.githubusercontent.com/wiki/analogdevicesinc/libiio/img/win_box.png)](https://ci.appveyor.com/api/projects/analogdevicesinc/libad9361-iio/artifacts/libad9361-setup.exe?branch=master) | Win32 : [![Latest 32-bit Windows zip](https://raw.githubusercontent.com/wiki/analogdevicesinc/libiio/img/win_box.png)](https://ci.appveyor.com/api/projects/analogdevicesinc/libad9361-iio/artifacts/libad9361-win32.zip?branch=master)<br />Win64: [![Latest 64-bit Windows zip](https://raw.githubusercontent.com/wiki/analogdevicesinc/libiio/img/win_box.png)](https://ci.appveyor.com/api/projects/analogdevicesinc/libad9361-iio/artifacts/libad9361-win64.zip?branch=master) |
| OS X                    |  [![OSX Status](https://api.travis-ci.org/analogdevicesinc/libad9361-iio.svg?branch=master&label=osx&passingTex=foo)](https://travis-ci.org/analogdevicesinc/libad9361-iio) |  OS X Mojave <br />(v 10.14) | [![OS-X package 10.14](https://raw.githubusercontent.com/wiki/analogdevicesinc/libiio/img/osx_box.png)](http://swdownloads.analog.com/cse/travis_builds/master_latest_libad9361-iio-darwin-10.14.4.pkg) | [![OS-X tarball 10.14](https://raw.githubusercontent.com/wiki/analogdevicesinc/libiio/img/osx_box.png)](http://swdownloads.analog.com/cse/travis_builds/master_latest_libad9361-iio-darwin-10.14.4.tar.gz) |
|                         |                     |  OS X High Sierra <br />(v 10.13) | [![OS-X package 10.13](https://raw.githubusercontent.com/wiki/analogdevicesinc/libiio/img/osx_box.png)](http://swdownloads.analog.com/cse/travis_builds/master_latest_libad9361-iio-darwin-10.13.6.pkg) | [![OS-X tarball 10.13](https://raw.githubusercontent.com/wiki/analogdevicesinc/libiio/img/osx_box.png)](http://swdownloads.analog.com/cse/travis_builds/master_latest_libad9361-iio-darwin-10.13.6.tar.gz) |
|                    |                     | macOS Sierra<br />(v 10.12) | [![OS-X package 10.12](https://raw.githubusercontent.com/wiki/analogdevicesinc/libiio/img/osx_box.png)](http://swdownloads.analog.com/cse/travis_builds/master_latest_libad9361-iio-darwin-10.12.6.pkg) | [![OS-X tarball 10.12](https://raw.githubusercontent.com/wiki/analogdevicesinc/libiio/img/osx_box.png)](http://swdownloads.analog.com/cse/travis_builds/master_latest_libad9361-iio-darwin-10.12.6.tar.gz) |
| Linux     | [![Linux Status](https://api.travis-ci.org/analogdevicesinc/libad9361-iio.svg?branch=master&label=linux)](https://travis-ci.org/analogdevicesinc/libad9361-iio) | Ubuntu Bionic Beaver<br />(v 18.04)<sup>1</sup>  | [![Debian](https://raw.githubusercontent.com/wiki/analogdevicesinc/libiio/img/deb.png)](http://swdownloads.analog.com/cse/travis_builds/master_latest_libad9361-iio-ubuntu-18.04-amd64.deb) | [![RPM File](https://raw.githubusercontent.com/wiki/analogdevicesinc/libiio/img/rpm.png)](http://swdownloads.analog.com/cse/travis_builds/master_latest_libad9361-iio-ubuntu-18.04-amd64.rpm) [![tar.gz](https://raw.githubusercontent.com/wiki/analogdevicesinc/libiio/img/linux_box.png)](http://swdownloads.analog.com/cse/travis_builds/master_latest_libad9361-iio-ubuntu-18.04-amd64.tar.gz) |
|  |  | Ubuntu Xenial Xerus<br />(v 16.04)<sup>1</sup> | [![Debian](https://raw.githubusercontent.com/wiki/analogdevicesinc/libiio/img/deb.png)](http://swdownloads.analog.com/cse/travis_builds/master_latest_libad9361-iio-ubuntu-16.04-amd64.deb) | [![RPM File](https://raw.githubusercontent.com/wiki/analogdevicesinc/libiio/img/rpm.png)](http://swdownloads.analog.com/cse/travis_builds/master_latest_libad9361-iio-ubuntu-16.04-amd64.rpm)  [![tar.gz file](https://raw.githubusercontent.com/wiki/analogdevicesinc/libiio/img/linux_box.png)](http://swdownloads.analog.com/cse/travis_builds/master_latest_libad9361-iio-ubuntu-16.04-amd64.tar.gz) |
|  |  | CentOS 7  | [![RPM File](https://raw.githubusercontent.com/wiki/analogdevicesinc/libiio/img/rpm.png)](http://swdownloads.analog.com/cse/travis_builds/master_latest_libad9361-iio-centos-7-x86_64.rpm)  | [![Debian](https://raw.githubusercontent.com/wiki/analogdevicesinc/libiio/img/deb.png)](http://swdownloads.analog.com/cse/travis_builds/master_latest_libad9361-iio-centos-7-x86_64.deb) [![tar.gz](https://raw.githubusercontent.com/wiki/analogdevicesinc/libiio/img/linux_box.png)](http://swdownloads.analog.com/cse/travis_builds/master_latest_libad9361-iio-centos-7-x86_64.tar.gz) |
|  |  | CentOS 6  | [![RPM File](https://raw.githubusercontent.com/wiki/analogdevicesinc/libiio/img/rpm.png)](http://swdownloads.analog.com/cse/travis_builds/master_latest_libad9361-iio-centos-6.10-x86_64.rpm) |  [![Debian](https://raw.githubusercontent.com/wiki/analogdevicesinc/libiio/img/deb.png)](http://swdownloads.analog.com/cse/travis_builds/master_latest_libad9361-iio-centos-6.10-x86_64.deb) [![tar.gz](https://raw.githubusercontent.com/wiki/analogdevicesinc/libiio/img/linux_box.png)](http://swdownloads.analog.com/cse/travis_builds/master_latest_libad9361-iio-centos-6.10-x86_64.tar.gz) |

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
