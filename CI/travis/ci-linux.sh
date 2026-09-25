#!/bin/bash

set -e

handle_centos() {
	local package=$1
	yum -y install yum-utils gcc
	yum config-manager --set-enabled powertools
	yum localinstall -y $package
	export CMAKE_OPTIONS="-DPYTHON_BINDINGS=ON -DENABLE_PACKAGING=ON .."
	export LD_LIBRARY_PATH=/usr/local/lib:/usr/local/lib64
}

handle_fedora() {
	local package=$1
	dnf install -y ./$package
	export LD_LIBRARY_PATH=/usr/local/lib:/usr/local/lib64
	export PKG_CONFIG_PATH=/usr/local/lib/pkgconfig:/usr/local/lib64/pkgconfig
	export CMAKE_OPTIONS="-DPYTHON_BINDINGS=ON -DENABLE_PACKAGING=ON .."
}

handle_default() {
	local package=$1
	DEBIAN_FRONTEND=noninteractive apt-get install -y rpm
	sudo dpkg -i $package
	export CMAKE_OPTIONS="-DPYTHON_BINDINGS=ON -DENABLE_PACKAGING=ON -DDEB_DETECT_DEPENDENCIES=ON .."
}

handle_opensuse() {
	local package=$1
	zypper in -y --allow-unsigned-rpm $package
	export CMAKE_OPTIONS="-DPYTHON_BINDINGS=ON -DENABLE_PACKAGING=ON .."
}

handle_"$1" "$2"

rm -f /usr/lib/python*/EXTERNALLY-MANAGED

# Build project
mkdir -p build
cd build
cmake $CMAKE_OPTIONS
sudo make && sudo make package && make test
sudo make install
ldconfig /usr/local/lib/ /usr/local/lib64/