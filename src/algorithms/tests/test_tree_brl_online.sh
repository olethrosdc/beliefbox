#! /bin/bash

ulimit -Sv 2000000
env=Chain
T=100
runs=1
for horizon in 1 2 3 4 5
do
	for state_samples in 1 2 3 4
	do
		for reward_samples in 1 2
		do
			for leaf in Min MeanMDP UBound
			do
				echo $horizon $state_samples $reward_samples $leaf 
				fname=tbrl-${env}-${horizon}h-${state_samples}s-${leaf}.out
				#sem -j 8
				./bin/online_algorithms --environment $env --n_episodes 1 --n_runs $runs --gamma 0.95 --n_steps $T --algorithm TBRL --horizon $horizon --max_samples $state_samples --max_reward_samples $reward_samples --leaf_value $leaf --reward_prior Beta &>$fname
			done
		done
	done
done

