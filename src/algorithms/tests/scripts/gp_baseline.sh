#! /bin/bash

n_states=8
runs=100
T=1000000
episodes=10000
steps=-1
gamma=0.999

for environment in Puddle #MountainCar Puddle
do

    outdir=$HOME/results/gaussian_processes/$environment/${n_states}
    mkdir -p $outdir

    for algorithm in QLearning Sarsa UCRL 
    do
        for eps in 01 #01 1
        do
            echo $algorithm $eps
            outfile=`mktemp`
            ./bin/online_algorithms --gamma $gamma --algorithm $algorithm --environment $environment --n_states $n_states --n_runs $runs --n_steps $T --episode_steps  $steps --n_episodes $episodes --epsilon 0.${eps} > $outfile
            grep EPISODE_RETURN $outfile >$outdir/${algorithm}_${eps}e.episode;
            grep REWARD $outfile >$outdir/${algorithm}_${eps}e.reward
            grep PAYOFF $outfile >$outdir/${algorithm}_${eps}e.payoff
        done
    done


    for s in 1 2  4 8
    do
        for algorithm in USampling LSampling
        do
            echo $algorithm $s
            outfile=`mktemp`
            ./bin/online_algorithms --gamma $gamma --algorithm $algorithm --environment $environment --n_states $n_states --n_runs $runs --n_steps $T --episode_steps $steps --n_episodes $episodes --epsilon 0.0 --max_samples $s  >$outfile
            grep EPISODE_RETURN $outfile >$outdir/${algorithm}_${s}s.episode;
            grep REWARD $outfile >$outdir/${algorithm}_${s}s.reward;
            grep PAYOFF $outfile >$outdir/${algorithm}_${s}s.payoff;
        done
    done
done
