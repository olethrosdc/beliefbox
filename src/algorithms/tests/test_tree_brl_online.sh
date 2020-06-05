#! /bin/bash

ulimit -Sv 2000000
env=Chain
for horizon in 0 1 2 3 
do
	for state_samples in 2 3 4 5
	do
		echo $horizon
		fname=tbrl-${env}-${horizon}h-${state_samples}s.out
		./bin/online_algorithms --environment $env --n_episodes 1 --n_runs 10 --gamma 0.99 --n_steps 10000 --algorithm TBRL --horizon $horizon --max_samples $state_samples >$fname&
	done
done

