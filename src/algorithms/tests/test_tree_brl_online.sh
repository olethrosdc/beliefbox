#! /bin/bash

ulimit -Sv 2000000
for horizon in 0 1 2 3 4 5 6 7 8
do
	for state_samples in 2 4 8
	do
		echo $horizon
		fname=tbrl-${horizon}h-${state_samples}s.out
		sem -j 4 ./bin/online_algorithms --environment DoubleLoop --n_states 5 --n_episodes 1 --n_runs 10 --gamma 0.99 --n_steps 10000 --algorithm TBRL --horizon $horizon --max_samples $state_samples >$fname
	done
done

