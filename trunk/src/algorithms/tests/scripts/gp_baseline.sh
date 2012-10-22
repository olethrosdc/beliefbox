#! /bin/bash

n_states=8
runs=1000
T=1000000
episodes=1000

for environment in MountainCar Puddle
do

    outdir=$HOME/results/gaussian_processes/$environment/${n_states}
    mkdir -p $outdir

    for algorithm in QLearning Sarsa UCRL 
    do
        for eps in 0 01 1
        do
            echo $algorithm $eps
            ./bin/online_algorithms --algorithm $algorithm --environment $environment --n_states $n_states --n_runs $runs --n_steps $T --n_episodes $episodes --epsilon 0.${eps} | grep EPISODE_RETURN >$outdir/${algorithm}_${eps}e.episode;
        done
    done


    for s in 1 2 4 8
    do
        for algorithm in USampling LSampling
        do
            echo $algorithm $s
            ./bin/online_algorithms --algorithm $algorithm --environment $environment --n_states $n_states --n_runs $runs --n_steps $T --n_episodes $episodes --epsilon 0.0 --max_samples $s | grep EPISODE_RETURN >$outdir/${algorithm}_${s}s.episode;
        done
    done
done
