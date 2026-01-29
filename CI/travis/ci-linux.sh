#!/bin/bash

set -x

handle_centos() {
	local package=$1
	yum -y install yum-utils gcc
	yum config-manager --set-enabled powertools
	yum localinstall -y $package
	export CMAKE_OPTIONS="-DPYTHON_BINDINGS=ON -DENABLE_PACKAGING=ON .."
}

handle_default() {
	local package=$1
	apt-get update
	DEBIAN_FRONTEND=noninteractive apt-get install -y build-essential devscripts debhelper dh-python python3 python3-setuptools
	sudo dpkg -i $package
	export USE_DEBUILD=true
}

handle_opensuse() {
	local package=$1
	zypper in -y --allow-unsigned-rpm $package
	export CMAKE_OPTIONS="-DPYTHON_BINDINGS=ON -DENABLE_PACKAGING=ON .."
}

handle_"$1" "$2"

python3 -m pip install pylibiio --no-binary :all:

if [ "$USE_DEBUILD" = "true" ]; then
	# Debian packaging with debuild
	# Extract version from CMakeLists.txt project() command
	VERSION=$(grep -oP 'project\s*\([^)]*VERSION\s+\K[0-9.]+' CMakeLists.txt)
	ARCHITECTURE=$(dpkg --print-architecture)
	SOURCE_DIR=$(grep "^Source:" packaging/debian/control | cut -d' ' -f2)
	CURRENT_DIR=$(basename "$PWD")

	# Replace placeholders in debian template files
	sed -i "s/@VERSION@/$VERSION-1/" packaging/debian/changelog
	sed -i "s/@DATE@/$(date -R)/" packaging/debian/changelog
	sed -i "s/@ARCHITECTURE@/$ARCHITECTURE/" packaging/debian/control

	# Create orig tarball
	pushd ..
	tar czf ${SOURCE_DIR}_${VERSION}.orig.tar.gz $CURRENT_DIR
	popd

	# Copy debian directory to root
	cp -r packaging/debian .
	# Ensure .install files are not executable (WSL filesystem issue)
	chmod -x debian/*.install

	# Build debian package
	debuild

	# Install the built packages
	sudo dpkg -i ../*.deb

	# Move packages to build directory
	mkdir -p build
	mv ../*.deb build/

else
	# CMake packaging for RPM distros
	mkdir -p build
	cd build
	cmake $CMAKE_OPTIONS
	sudo make && sudo make package && make test
	sudo make install
	cd ..
fi

ldconfig

cd bindings/python
pip install -r requirements_dev.txt
export LD_LIBRARY_PATH=/usr/lib/:/usr/lib64/:/usr/local/lib/:/usr/local/lib64/
python3 -m pip install pytest
python3 -m pytest -vs --skip-scan
