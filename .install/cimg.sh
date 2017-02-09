#!/bin/bash
# cimg.sh
#   Script to install CImg-1.6.1

## Parse command line arguments to get install PREFIX
SHELL_DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
source $SHELL_DIR/parse.sh

if [ -d "$PREFIX/include" ]; then
  echo "Detected CImg installation in $PREFIX"
  exit 0
fi

which wget || ( echo Cannot auto-install cimg without wget && exit 1 )
echo making $PREFIX/include
mkdir -p $PREFIX/include
echo Downloading CImg
wget http://cimg.eu/files/CImg_1.6.4.zip || exit 1
unzip CImg_1.6.4.zip || exit 1
mv CImg-1.6.4/CImg.h ${PREFIX}/include/CImg.h || exit 1
rm ./CImg_1.6.4.zip || exit 1
rm -rf CImg-1.6.4 || exit 1
