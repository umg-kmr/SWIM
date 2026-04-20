#!/bin/bash

# Exit the script if any command fails
set -e

echo "Starting build..."

#Change these paths to the location of boost libraries and SWIM package
export BOOST_PATH=/mnt/misc/boost_1_90_0 #/path/to/boost/libraries
export SWIM_PATH=/mnt/misc/SWIM_v1 #/path/to/SWIM

#prints the set library paths
echo "BOOST path set to: $BOOST_PATH"
echo "SWIM path set to: $SWIM_PATH"

# Compile C++ files
# You can adjust the compilation flags as needed

echo "compiling background module..."
g++ -shared -fPIC -I $BOOST_PATH -O3 -march=native -mtune=native -ftree-vectorize -funroll-loops -o $SWIM_PATH/GQ_Calculator/bg/libbg.so $SWIM_PATH/GQ_Calculator/bg/model_calc.cpp -lm -fopenmp

echo "compiling perturbation module..."
g++ -shared -fPIC -I $BOOST_PATH -O3 -march=native -mtune=native -ftree-vectorize -funroll-loops -o $SWIM_PATH/GQ_Calculator/pert/libpert.so $SWIM_PATH/GQ_Calculator/pert/model_calc.cpp -lm -fopenmp

echo "compiling power spectrum module..."
g++ -shared -fPIC -I $BOOST_PATH -O3 -march=native -mtune=native -ftree-vectorize -funroll-loops -o $SWIM_PATH/PS_Calculator/libmodel.so $SWIM_PATH/PS_Calculator/model_calc.cpp -lm -fopenmp

echo "compiling semi-analytical power spectrum module..."
g++ -shared -fPIC -I $BOOST_PATH -O3 -march=native -mtune=native -ftree-vectorize -funroll-loops -o $SWIM_PATH/SA_PS_Calculator/libmodel.so $SWIM_PATH/SA_PS_Calculator/model_calc.cpp -lm -fopenmp

echo "compilation finished."

