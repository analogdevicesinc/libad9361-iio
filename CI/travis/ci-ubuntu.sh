#!/bin/bash

set -ex
uname -a
DEBIAN_FRONTEND=noninteractive apt install -y make graphviz libaio-dev \
	libavahi-client-dev libavahi-common-dev libusb-1.0-0-dev \
	rpm tar bzip2 gzip libserialport-dev python3-pip
dpkg -i /ci/build/*.deb
# refresh the loader cache so the just-installed libiio is discoverable
# (its CPack .deb doesn't run ldconfig) before pylibiio's source build looks for it
ldconfig

# pip on PEP 668 distros (Debian 12+/Ubuntu 24.04+) needs --break-system-packages;
# older pip (Ubuntu 22.04) doesn't know the flag, so fall back to a plain install.
pip_install() { python3 -m pip install --break-system-packages "$@" 2>/dev/null || python3 -m pip install "$@"; }

pip_install setuptools wheel
pip_install pylibiio
pip_install sphinx
pip_install sphinx-rtd-theme furo

echo "$PWD"

mkdir -p build
cd build
cmake -DPYTHON_BINDINGS=ON -DENABLE_PACKAGING=ON -DDEB_DETECT_DEPENDENCIES=ON -DWITH_DOC=OFF \
	-DCMAKE_BUILD_TYPE="${CMAKE_BUILD_TYPE:-Release}" ..
make && make package && make test
make install
ldconfig
export LD_LIBRARY_PATH=/usr/local/lib:${LD_LIBRARY_PATH:-}
cd ../bindings/python
pip_install -r requirements_dev.txt
python3 -m pytest -vs --skip-scan
