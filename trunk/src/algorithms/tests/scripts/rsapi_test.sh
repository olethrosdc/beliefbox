#! /bin/bash

resdir=$HOME/results/rsapi
gamma=0.95
n_iter=25
horizon=1000
for n_states in 10 20 50 100 200 500 1000
do
    for n_rollouts in 10 100 1000
    do
        for environment in MountainCar Pendulum
        do
            outdir=$resdir/$environment/
            mkdir -p $outdir
            dst=$outdir/s${n_states}_r${n_rollouts}
            cmd="$HOME/projects/beliefbox/src/algorithms/tests/bin/rsapi --environment $environment --n_states $n_states --n_rollouts $n_rollouts --gamma $gamma --horizon $horizon --n_iter $n_iter"
            echo ${cmd} >$dst.cmd
            for run in 1 2 3 4 5 6 7 8 
            do
                nice -n 19 ${cmd} >${dst}_run${run}.out &
            done
            wait;
        done
    done
done