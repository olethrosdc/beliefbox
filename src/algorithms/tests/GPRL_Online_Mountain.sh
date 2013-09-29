#! /bin/bash

make DBG_OPT=DBG -j8  GPRL

for((n_runs = 1; n_runs<=100; n_runs = n_runs + 1))
do
	outdir=./GPRLMountainCarOnline/$n_rollout/$n_runs
	mkdir -p $outdir
	./bin/GPRL --online --environment MountainCar --discount 0.999 --n_training_episodes 1000 --n_training_e_steps 1000 --grid 4 --scale 1.0 --seed_file ????
	mv online_GPRL_RESULTS_STEPS_MountainCar $outdir
	mv online_GPRL_RESULTS_REWARDS_MountainCar $outdir
done

for((n_runs = 1; n_runs<=100; n_runs = n_runs + 1))
do
    outdir=./GPRLMountainCarOnline/$n_rollout/$n_runs
    mkdir -p $outdir
    ./bin/GPRL --online --environment MountainCar --discount 0.999 --n_training_episodes 1000 --n_training_e_steps 1000 --grid 4 --scale 1.0 --seed_file ????
    mv online_GPRL_RESULTS_STEPS_MountainCar $outdir
    mv online_GPRL_RESULTS_REWARDS_MountainCar $outdir
done


