#!/bin/sh
#
# Installs the lemon library in the current directory.

# display commands as they are executed
set -x

# remove potentially existing lemon files
rm -rf lib
rm -rf include
rm -rf bin

# extract lemon archive
tar -xzf lemon-1.2.3.tar.gz

# save current directory for installation prefix
lemon_prefix=`pwd`

# build
cd lemon-1.2.3
./configure --prefix=$lemon_prefix
make
make install

# remove source
cd ..
rm -rf lemon-1.2.3
