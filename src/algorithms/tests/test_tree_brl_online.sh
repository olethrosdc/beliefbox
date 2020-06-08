#! /bin/bash

ulimit -Sv 2000000
env=Chain
T=1000
for horizon in 1 2 3 4 5
do
	for state_samples in 1 2 3 
	do
		for leaf in Min Mean UBound LBound
		do
			echo $horizon $state_samples
			fname=tbrl-${env}-${horizon}h-${state_samples}s-${leaf}.out
			sem -j 8 ./bin/online_algorithms --environment $env --n_episodes 1 --n_runs 10 --gamma 0.95 --n_steps $T --algorithm TBRL --horizon $horizon --max_samples $state_samples --leaf_value $leaf &>$fname
		done
	done
done

