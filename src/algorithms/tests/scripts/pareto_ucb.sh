#! /bin/bash

exc=$HOME/projects/beliefbox/src/algorithms/tests/bin/pareto_ucb
actions=16
outcomes=4
horizon=1000000
runs=10

for actions in 2 4 8 16 32 64 128
do
	resdir=$HOME/results/pareto_ucb/${actions}A_${outcomes}S_${horizon}T_WUCB/
	mkdir -p $resdir
	results=$resdir/regret.out
	errors=$resdir/errors.out
	echo $actions $outcomes $horizon $runs $epsilon
	time ${exc} $actions $outcomes $horizon $runs 2 0.0 > $results 2> $errors

	resdir=$HOME/results/pareto_ucb/${actions}A_${outcomes}S_${horizon}T_HUCB/
	mkdir -p $resdir
	results=$resdir/regret.out
	errors=$resdir/errors.out
	echo $actions $outcomes $horizon $runs $epsilon
	time ${exc} $actions $outcomes $horizon $runs 3 0.0 > $results 2> $errors
	
	for epsilon in 0.0 0.1 1.0 #0.01 0.1 0.5 1.0
	do
		resdir=$HOME/results/pareto_ucb/${actions}A_${outcomes}S_${horizon}T_${epsilon}E/
		mkdir -p $resdir
		results=$resdir/regret.out
		errors=$resdir/errors.out
		echo $actions $outcomes $horizon $runs $epsilon
		time ${exc} $actions $outcomes $horizon $runs 1 $epsilon > $results 2> $errors
	done
done