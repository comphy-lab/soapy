#!/bin/bash

# Save the original directory
ORIG_DIR=$(pwd)

mkdir -p $1

cp $1.c $1/
cd $1


CC99='mpicc -std=c99' qcc -Wall -O2 -D_MPI=1 -disable-dimensions -I${ORIG_DIR}/src-local -I${ORIG_DIR}/../src-local $1.c -o $1 -lm

mpirun -np $2 ./$1