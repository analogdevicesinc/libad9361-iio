#!/bin/bash

set -x
uname -a

ARCH=$(uname -m)

# For arm32 only: save the pre-installed cmake (3.23+) before apt overwrites it with 3.16
if [ "$ARCH" = "armv7l" ]; then
    ORIGINAL_CMAKE=$(which cmake 2>/dev/null)
    echo "arm32 detected - Original cmake: $ORIGINAL_CMAKE"
fi

apt-get update
DEBIAN_FRONTEND=noninteractive apt-get install -y build-essential make devscripts debhelper \
	dh-python python3 python3-setuptools python3-pip \
	libaio-dev libavahi-client-dev libavahi-common-dev libusb-1.0-0-dev \
	libserialport-dev graphviz cmake \
	rpm tar bzip2 gzip

# For arm32 only: restore the original cmake (3.23+) instead of apt's 3.16
# CMake 3.16 from apt has a bug with compiler detection under QEMU arm32
if [ "$ARCH" = "armv7l" ] && [ -n "$ORIGINAL_CMAKE" ] && [ -x "$ORIGINAL_CMAKE" ]; then
    ln -sf "$ORIGINAL_CMAKE" /usr/bin/cmake
    echo "Restored cmake from: $ORIGINAL_CMAKE"
fi

# Install libiio from pre-built packages
dpkg -i /ci/build/*.deb || true
# Configure any pending packages and fix broken dependencies
dpkg --configure -a
apt-get install -f -y
ldconfig

python3 -m pip install pylibiio --no-binary :all:
python3 -m pip install sphinx
python3 -m pip install sphinx-rtd-theme

echo "$PWD"

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

# Copy debian directory to root
cp -r packaging/debian .
# Ensure .install files are not executable (WSL filesystem issue)
chmod -x debian/*.install

# Create orig tarball
pushd ..
tar czf ${SOURCE_DIR}_${VERSION}.orig.tar.gz $CURRENT_DIR
popd

# Build debian package
debuild

# Install the built packages
sudo dpkg -i ../*.deb

# Move packages to build directory
mkdir -p build
mv ../*.deb build/

ldconfig

cd bindings/python

pip3 install -r requirements_dev.txt

python3 -m pytest -vs --skip-scan
