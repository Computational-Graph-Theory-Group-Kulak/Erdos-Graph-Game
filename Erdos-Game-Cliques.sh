#!/bin/bash

dir=/home/u0149846/algorithms/erdos-stijn/$1-vertices
edges=$((($1*($1-1)+3)/4))

echo "edges: $edges"

mkdir -p $dir

/home/u0149846/nauty2_8_6/geng $1 $edges 2> $dir/nauty.e$i | /home/u0149846/algorithms/Erdos-Stijn/./ErdosGame $1 $(($1*($1-1)/2)) 2> $dir/error.txt 1>$dir/results.txt
echo "Done"