#!/bin/bash
echo ./$1-vertices/$2/$3/$4
dir=./$1-vertices/$2/$3/$4

mkdir -p $dir
./gen-colex $1 $2 $3 $4 -b$5 -p$7 2> $dir/gen_error.txt | ./erdos-solver $1 -m$6 $3 $4 -g -b$5 -p$7 2> $dir/solv_error.txt 1>$dir/results.txt
echo "Done"