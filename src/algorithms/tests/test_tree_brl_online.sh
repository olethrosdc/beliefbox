#! /bin/bash

ulimit -Sv 2000000
env=Chain
T=1000
for horizon in 0 1 2 3 4 5
do
	for state_samples in 1 3 #2 
	do
		echo $horizon $state_samples
		fname=tbrl-${env}-${horizon}h-${state_samples}s.out
		sem -j 8 ./bin/online_algorithms --environment $env --n_episodes 1 --n_runs 10 --gamma 0.95 --n_steps $T --algorithm TBRL --horizon $horizon --max_samples $state_samples &>$fname
	done
done

