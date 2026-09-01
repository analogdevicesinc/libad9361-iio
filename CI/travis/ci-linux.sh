#!/bin/bash

set -x

handle_centos() {
	local package=$1
	yum -y install yum-utils gcc
	yum config-manager --set-enabled powertools
	yum localinstall -y $package
	export CMAKE_OPTIONS="-DPYTHON_BINDINGS=ON -DENABLE_PACKAGING=ON .."
	export LD_LIBRARY_PATH=/usr/local/lib64/
}

handle_default() {
	local package=$1
	DEBIAN_FRONTEND=noninteractive apt-get install -y rpm
	dpkg -i $package
	export CMAKE_OPTIONS="-DPYTHON_BINDINGS=ON -DENABLE_PACKAGING=ON -DDEB_DETECT_DEPENDENCIES=ON .."
}

handle_opensuse() {
	local package=$1
	zypper in -y --allow-unsigned-rpm $package
	export CMAKE_OPTIONS="-DPYTHON_BINDINGS=ON -DENABLE_PACKAGING=ON .."
}

handle_"$1" "$2"

python3 -m pip install pylibiio --no-binary :all:
# Build project
mkdir -p build
cd build
cmake $CMAKE_OPTIONS
make && make package && make test
make install
ldconfig
cd ..
cd bindings/python
uv pip install -r requirements_dev.txt
uv pip install 'pytest==7.2.0'
python3 -m pytest -vs --skip-scan
