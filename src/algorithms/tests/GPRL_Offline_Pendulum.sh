#! /bin/bash

##make DBG_OPT=DBG -j8  GPRL

n_cpu=48

resdir=$HOME/results/GPRL
for((n_rollout=10;n_rollout<=45;n_rollout = n_rollout+5)) 
do
    for((n_runs = 1; n_runs<=100; n_runs = n_runs + 1))
    do
	outdir=$resdir/GPRLPendulumOffline/$n_rollout/$n_runs
	mkdir -p $outdir
	cd $outdir
	sem -j 8 --id gprl_mutex GPRL --environment Pendulum --discount 0.95 --n_training_episodes $n_rollout --n_training_e_steps 40 --n_testing_episodes 1000 --n_testing_e_steps 3000 --grid 3 --scale 0.5 --seed $n_runs --seed_file $HOME/data/r1e7.bin 2>&1 >out.log
    done
done

for((n_rollout=50;n_rollout<=1000;n_rollout = n_rollout+50)) 
do
    for((n_runs = 1; n_runs<=100; n_runs = n_runs + 1))
    do
	outdir=$resdir/GPRLPendulumOffline/$n_rollout/$n_runs
	mkdir -p $outdir
	cd $outdir
	sem -j 8 --id gprl_mutex GPRL --environment Pendulum --discount 0.95 --n_training_episodes $n_rollout --n_training_e_steps 40 --n_testing_episodes 1000 --n_testing_e_steps 3000 --grid 3 --scale 0.5 --seed $n_runs --seed_file $HOME/data/r1e7.bin 2>&1 >out.log
    done
done

