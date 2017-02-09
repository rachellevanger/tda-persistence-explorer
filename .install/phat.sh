#!/bin/bash
# phat.sh
#   Script to install PHAT

## Parse command line arguments to get install PREFIX
SHELL_DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
source $SHELL_DIR/parse.sh

if [ -d "$PREFIX/include" ]; then
  echo "Detected phat installation in $PREFIX"
  exit 0
fi

echo Downloading PHAT
git clone https://github.com/blazs/phat.git || exit 1

echo making $PREFIX/include
mkdir -p $PREFIX/include || exit 1

echo Installing PHAT headers to $PREFIX/include 
cp -rf phat/include/phat $PREFIX/include/phat || exit 1

echo Removing phat 
rm -rf phat || exit 1
