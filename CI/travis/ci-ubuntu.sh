#!/bin/bash

set -ex
uname -a
DEBIAN_FRONTEND=noninteractive apt install -y make graphviz libaio-dev \
	libavahi-client-dev libavahi-common-dev libusb-1.0-0-dev \
	rpm tar bzip2 gzip libserialport-dev
dpkg -i /ci/build/*.deb
# refresh the loader cache so the just-installed libiio is discoverable
# (its CPack .deb doesn't run ldconfig) before pylibiio's source build looks for it
ldconfig

# The caller activates a venv (PEP 668 distros forbid system-wide pip installs).
pip install setuptools wheel pylibiio sphinx sphinx-rtd-theme furo

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
pip install -r requirements_dev.txt
python3 -m pytest -vs --skip-scan
