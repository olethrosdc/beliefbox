#! /bin/bash

n_runs=1000
n_episodes=10000
n_steps=10000
gamma=0.99
lambda=0.8
n_actions=2
for environment in MountainCar Pendulum #MountainCar Pendulum
do
    for random in 0 0.01
    do
        for K in 2 3 4 5
        do
            algorithm=QLearning
            resdir=$HOME/results/$environment/$algorithm/rand${random}
            mkdir -p $resdir
            echo ./bin/pomdp_algorithms $K $n_actions $gamma $lambda $random $n_runs $n_episodes $n_steps $algorithm $environment | tee $resdir/run.params
            nice -n 19 ./bin/pomdp_algorithms $K $n_actions $gamma $lambda $random $n_runs $n_episodes $n_steps $algorithm $environment | grep REWARD > ${resdir}/K${K}reward.out &
            algorithm=BVMM
            resdir=$HOME/results/$environment/$algorithm/rand${random}
            mkdir -p $resdir
            for D in 0 1 2 3 4
            do
                echo ./bin/pomdp_algorithms $K $n_actions $gamma $lambda $random $n_runs $n_episodes $n_steps $algorithm $environment $D | tee $resdir/run.params
                nice -n 19 ./bin/pomdp_algorithms $K $n_actions $gamma $lambda $random $n_runs $n_episodes $n_steps $algorithm $environment $D | grep REWARD > ${resdir}/K${K}_D${D}_reward.out &
            done
            wait;
        done
    done
done
