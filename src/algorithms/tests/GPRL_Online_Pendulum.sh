#! /bin/bash

make DBG_OPT=DBG -j8  GPRL

for((n_runs = 1; n_runs<=100; n_runs = n_runs + 1))
do
    outdir=./GPRLPendulumOnline/$n_rollout/$n_runs
    mkdir -p $outdir
    ./bin/GPRL --online --environment Pendulum --discount 0.95 --n_training_episodes 1000 --n_training_e_steps 3000 --grid 3 --scale 0.5 --seed_file ????
    mv online_GPRL_RESULTS_STEPS_Pendulum $outdir
    mv online_GPRL_RESULTS_REWARDS_Pendulum $outdir
done

for((n_runs = 1; n_runs<=100; n_runs = n_runs + 1))
do
    outdir=./GPRLPendulumOnline/$n_rollout/$n_runs
    mkdir -p $outdir
    ./bin/GPRL --online --environment Pendulum --discount 0.95 --n_training_episodes 1000 --n_training_e_steps 3000 --grid 3 --scale 0.5 --seed_file ????
    mv online_GPRL_RESULTS_STEPS_Pendulum $outdir
    mv online_GPRL_RESULTS_REWARDS_Pendulum $outdir
done


