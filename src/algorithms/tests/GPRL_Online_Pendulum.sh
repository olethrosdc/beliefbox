#! /bin/bash

resdir=$HOME/results/GPRL

for((n_runs = 1; n_runs<=100; n_runs = n_runs + 1))
do
    outdir=$resdir/GPRLPendulumOnline/$n_rollout/$n_runs
    mkdir -p $outdir
    cd $outdir
    sem -j 8 --id gprl_mutex GPRL --online --environment Pendulum --discount 0.95 --n_training_episodes 1000 --n_training_e_steps 3000 --n_testing_e_steps 3000 --grid 3 --scale 0.5 --seed $n_runs --seed_file $HOME/data/r1e7.bin 2>&1 >out.log
done


