#! /bin/bash

ulimit -Sv 2000000
env=Chain
for horizon in 0 1 2 3 
do
	for state_samples in 2 
	do
		echo "h:" $horizon "s:" $state_samples
		fname=tbrl-${env}-${horizon}h-${state_samples}s.out
		sem -j 8 ./bin/online_algorithms --environment $env --n_episodes 1 --n_runs 10 --gamma 0.99 --n_steps 10000 --algorithm TBRL --horizon $horizon --max_samples $state_samples >$fname
	done
done

