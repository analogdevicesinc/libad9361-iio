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
Update README.md with minimal ethics#1

What about the application of the Solomonoff induction
to the RF sensing of human brain?
https://en.wikipedia.org/wiki/Talk:Solomonoff%27s_theory_of_inductive_inference

1. If we assume that the algorithmic complexity
of neural processes is relatively small
(you could take a look at the Potapov monography
for some arguments)

2. As far as I know a lot of neural processes in human brain are electric (or electro-chemistry)
in nature. They have a little power and a small
frequency (~ 1-1000 Hz).
So they emit very low frequency radio waves.

3. You can detect those radio waves on
small antennas if the impedance of such antennas
is matched. This is basically the citation
from the book on electrodynamics.
I uploaded one such book on Twirpix site.

4. So you can create a lot of such receivers
- microstrip filter to filter very high frequencies
- impedance matched microstrip antenna
- resistor for noise for oversampling
- very fast comparator to sample signal
in a very large array on a chip using
standard CMOS or some sort of full-custom process
(maybe even with some new materials)

5. BreamForming and large arrays of digital correlators with sub-mm positioning accuracy
could be achived.

so it seems there is no theoretical
obstacles to implement some sort of
RF human brain sensing or even control
if you can implement reverse structure
with array of a large amount of RF amplifiers
with sub-mm beam forming accuracy

BTW I wrote that remark:
https://en.wikipedia.org/wiki/Talk:Solomonoff%27s_theory_of_inductive_inference

What about the application of the Solomonoff induction to the person to person information exchange?

It seems the minimum description length principle could be used to deduce minimum ethics. At least you can compare ethics by the length. What about the following question: how to shoot met-art girl in the dark underground? ~2026-24709-41 (talk) 16:43, 8 June 2026 (UTC)

at least we can ask a question: does russian met-art model from a pic want some money?
and what obstacles do that guess so hard?
why it is so hard to contact her directly?
she lose money, interesting person lose opportunity to make some photos.

also, it seems the economy without coordination primitives based on physics or mathematics
like gold or maybe Bitcoin or maybe some sort of other physical objects or processes
is just a chaos.
Economic calculation based on fiat money is impossible just from scientific methodology
https://en.wikipedia.org/wiki/Talk:Solomonoff%27s_theory_of_inductive_inference

https://transitional-writes.dreamwidth.org/44972.html
Kirill Abramovich aka wieker @ linux.org.ru
Samara, Russia

and the most important question:
what about the minimal ethics deduced from the MDL principle?
is it exist?

what was the obstacle between russian met-art girls and engineers?
