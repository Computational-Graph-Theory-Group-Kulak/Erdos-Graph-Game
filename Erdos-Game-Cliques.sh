#!/bin/bash

dir=./$1-vertices

edges=$((($1*($1-1)+3)/4))

echo "edges: $edges"

mkdir -p $dir

./nauty2_8_6/geng $1 $edges 2> $dir/error_geng.txt 1> $dir/endconfigurations.e

counts=$(wc -l < $dir/endconfigurations.e)
echo "counts: $counts"

echo "$counts
$(cat $dir/endconfigurations.e)" > $dir/endconfigurations.e

cat $dir/endconfigurations.e | ./erdos-solver $1 2> $dir/error_solv.txt 1>$dir/results.txt
echo "Done"