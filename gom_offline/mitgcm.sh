#!/bin/bash

# ---------- #
# Creates directories, files, and executable to run MITgcm
# Written by: Taylor Shropshire
# Date: 6/6/16
# --------- #

# Remove and Create Working Directories
/bin/rm -rf ./build ./run
/bin/mkdir ./build ./run

# Compile MITgcm
cd ./build
echo "Running genmake..."
sleep 1
../../../tools/genmake2 -mods=../code
echo "Make clean"
sleep 1
make clean
echo "Make depend"
sleep 1
make depend
echo "Final make"
sleep 1
make


# Organize input files and the MITgcm exe. created above
cd ../run
ln -s ../input/* .
ln -s ../build/mitgcmuv

echo "Now switch to <./run> directory"
