#!/bin/bash

set -x

yum -y install yum-utils gcc
yum config-manager --set-enabled powertools

cd /opt
wget https://www.python.org/ftp/python/3.5.1/Python-3.5.1.tgz
tar xzf Python-3.5.1.tgz
cd Python-3.5.1
./configure
make altinstall
cd /ci

yum -y install epel-release bzip2 gzip

yum localinstall -y /ci/build/*.rpm

python3.5 -m pip install pylibiio --no-binary :all:

# Build project
sudo ln -fs /usr/bin/python3.5 /usr/bin/python
python3.5 -m ensurepip
mkdir -p build
cd build
cmake -DPYTHON_BINDINGS=ON -DENABLE_PACKAGING=ON ..
sudo make && sudo make package && make test
sudo make install
ldconfig
cd ..
cd bindings/python
pip install -r requirements_dev.txt
export LD_LIBRARY_PATH=/usr/local/lib/
python3.5 -m pip install pytest
python3.5 -m pytest -vs --skip-scan
