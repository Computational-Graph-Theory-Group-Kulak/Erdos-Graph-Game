#!/bin/bash
dir=./$1-vertices/Compare

edges=$((($1*($1-1)+3)/4))

mkdir -p $dir
./gen-colex $1 2> $dir/error_generator.txt 1> $dir/results_generator.txt
./nauty2_8_6/geng $1 $edges 2> $dir/error_geng.txt 1> $dir/results_geng.txt
counts_geng=$(wc -l < $dir/results_geng.txt)
read -r counts_generator<$dir/results_generator.txt

echo "geng generated $counts_geng non-isomorphic end configurations and the generator created $counts_generator"