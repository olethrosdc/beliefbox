#! /bin/bash

n_runs=10
n_episodes=100
n_steps=1000
gamma=0.99
lambda=0.8
n_actions=2
for algorithm in BVMM #QLearning #Sarsa
do
    for environment in MountainCar Pendulum
    do
        for random in 0 0.001 0.01 0.1
        do
            resdir=$HOME/results/$environment/$algorithm/rand${random}
            mkdir -p $resdir
            for K in 2 3 4 5
            do
                echo ./bin/pomdp_algorithms $K $n_actions $gamma $lambda $random $n_runs $n_episodes $n_steps $algorithm $environment | tee $resdir/run.params
                nice -n 19 ./bin/pomdp_algorithms $K $n_actions $gamma $lambda $random $n_runs $n_episodes $n_steps $algorithm $environment | grep REWARD > ${resdir}/K${K}reward.out 2>/dev/null;
            done
            wait;
            for K in 6 7 8 10
            do
                echo ./bin/pomdp_algorithms $K $n_actions $gamma $lambda $random $n_runs $n_episodes $n_steps $algorithm $environment | tee $resdir/run.params
                nice -n 19 ./bin/pomdp_algorithms $K $n_actions $gamma $lambda $random $n_runs $n_episodes $n_steps $algorithm $environment | grep REWARD > ${resdir}/K${K}reward.out 2>/dev/null;
            done
            wait;
        done
    done
done
