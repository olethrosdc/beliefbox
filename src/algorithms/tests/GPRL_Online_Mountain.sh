#! /bin/bash

for((n_runs = 1; n_runs<=100; n_runs = n_runs + 1))
do
	outdir=./GPRLMountainCarOnline/$n_rollout/$n_runs
	mkdir -p $outdir
	./bin/GPRL --online --environment MountainCar --discount 0.999 --n_training_episodes 1000 --n_training_e_steps 1000 --grid 4 --scale 1.0 --seed_file $HOME/data/r1e7.bin
done

for((n_runs = 1; n_runs<=100; n_runs = n_runs + 1))
do
    outdir=./GPRLMountainCarOnline/$n_rollout/$n_runs
    mkdir -p $outdir
    sem -j 8 --id gprl_mutex GPRL --online --environment MountainCar --discount 0.999 --n_training_episodes 1000 --n_training_e_steps 1000 --grid 4 --scale 1.0 --seed_file $HOME/data/r1e7.bin  2>&1 >out.log
done


