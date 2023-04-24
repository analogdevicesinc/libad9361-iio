#!/bin/bash

set -x
uname -a
DEBIAN_FRONTEND=noninteractive apt install -y make graphviz libaio-dev \
	libavahi-client-dev libavahi-common-dev libusb-1.0-0-dev \
	rpm tar bzip2 gzip libserialport-dev python3-pip
dpkg -i /ci/build/*.deb
python3 -m pip install pylibiio --no-binary :all:
python3 -m pip install sphinx
python3 -m pip install sphinx-rtd-theme furo

echo "$PWD"

mkdir -p build
cd build
cmake -DPYTHON_BINDINGS=ON -DENABLE_PACKAGING=ON -DDEB_DETECT_DEPENDENCIES=ON -DWITH_DOC=OFF ..
make && make package && make test
make install
ldconfig
cd ../bindings/python
pip3 install -r requirements_dev.txt
python3 -m pytest -vs --skip-scan
