#! /bin/bash

summary=test_tree_brl.dat
rm -f $summary
for horizon in 0 1 2 4 8
do
	for leaf in 0 2 3 4 5
	do
		fname=out-${horizon}-${leaf}
		./bin/test_tree_brl $horizon $leaf >$fname
		grep H $fname >> $summary
	done
done
	
