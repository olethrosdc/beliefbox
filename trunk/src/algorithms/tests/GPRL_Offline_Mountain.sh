#! /bin/bash

make DBG_OPT=DBG -j8  GPRL

for((n_rollout=10;n_rollout<=45;n_rollout = n_rollout+5)) 
do
    for((n_runs = 1; n_runs<=100; n_runs = n_runs + 1))
    do
	outdir=./GPRLMountainCarOffline/$n_rollout/$n_runs
	mkdir -p $outdir
	./bin/GPRL --environment MountainCar --discount 0.999 --n_training_episodes $n_rollout --n_training_e_steps 40 --n_testing_episodes 1000 --n_testing_e_steps 1000 --grid 4 --scale 1.0 --seed_file ????
	mv offline_GPRL_RESULTS_STEPS_MountainCar $outdir
	mv offline_GPRL_RESULTS_REWARDS_MountainCar $outdir
    done
done

for((n_rollout=50;n_rollout<=1000;n_rollout = n_rollout+50)) 
do
    for((n_runs = 1; n_runs<=100; n_runs = n_runs + 1))
    do
	outdir=./GPRLMountainCarOffline/$n_rollout/$n_runs
	mkdir -p $outdir
	./bin/GPRL --environment MountainCar --discount 0.999 --n_training_episodes $n_rollout --n_training_e_steps 40 --n_testing_episodes 1000 --n_testing_e_steps 1000 --grid 4 --scale 1.0 --seed_file ????
	mv offline_GPRL_RESULTS_STEPS_MountainCar $outdir
	mv offline_GPRL_RESULTS_REWARDS_MountainCar $outdir
    done
done

