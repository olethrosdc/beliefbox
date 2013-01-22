#! /bin/bash

resdir=$HOME/results/rsapi
gamma=0.999
n_iter=25
horizon=1000
for n_states in 10 20 50 100 200 500 1000
do
    for n_rollouts in 10 20 50 100 200 #1000
    do
        for environment in Pendulum MountainCar
        do
            for method in uniform error_bound upper_bound
	    do
            echo "$environment s:$n_states n:$n_rollouts"
            outdir=$resdir/$environment/${method}
            mkdir -p $outdir
            dst=$outdir/s${n_states}_r${n_rollouts}
            cmd="$HOME/projects/beliefbox/src/algorithms/tests/bin/rsapi --environment $environment --n_states $n_states --n_rollouts $n_rollouts --gamma $gamma --horizon $horizon --n_iter $n_iter --${method}"
            echo ${cmd} >$dst.cmd
            nice -n 19 ${cmd} >${dst}_run${run}.out &
            nice -n 19 ${cmd} --group >${dst}_run${run}_group.out &
            nice -n 19 ${cmd} --group --resample >${dst}_run${run}_group_resample.out &

            nice -n 19 ${cmd} --lipschitz 1 >${dst}_run${run}_lipschitz.out &
            nice -n 19 ${cmd} --group --lipschitz 1 >${dst}_run${run}_group_lipschitz.out &
            nice -n 19 ${cmd} --group --resample --lipschitz 1 >${dst}_run${run}_group_resample_lipschitz.out &

            wait;
        done
done
    done
done
