#! /bin/bash

ulimit -Sv 4000000
n=100
for horizon in 0 1 2 3 3 4 5 6 8 
do
	for leaf in 0 2 3 4 5
	do
		fname=out-${horizon}-${leaf}
		sem -j6 ./bin/test_tree_brl ${horizon} ${leaf} $n >$fname 
	done
done

