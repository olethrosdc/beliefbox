#! /bin/bash

ulimit -Sv 4000000
n=10
T=100000
gamma=0.95
for horizon in 0 1 2 3 3 4 5 6 8 
do
	for leaf in 0 2 3 4 5
	do
		fname=tbrl-${horizon}-${leaf}
		sem -j8 ./bin/test_tree_brl ${horizon} ${leaf} $n $T $gamma >$fname.out
#		grep "# reward" $fname.out >$fname.reward
	done
done

