#!/bin/bash

set -x

DEBIAN_FRONTEND=noninteractive apt-get install -y rpm

sudo dpkg -i /ci/build/*.deb

python3 -m pip install pylibiio --no-binary :all:

# Build project
mkdir -p build
cd build
sudo cmake -DPYTHON_BINDINGS=ON -DENABLE_PACKAGING=ON -DDEB_DETECT_DEPENDENCIES=ON ..
sudo make && make package && make test
sudo make install
ldconfig
cd ..
cd bindings/python
pip3 install -r requirements_dev.txt
python3 -m pytest -vs --skip-scan
