#! /bin/bash

n_actions=2
gamma=0.99
n_runs=1000
n_episodes=10000
n_steps=1000
exec_dir=~/projects/beliefbox/src/algorithms/tests/bin/

for n_states in 4 8 16 32 64 128
do
    res_dir=~/experiments/prior/dirichlet/${n_states}s_${n_actions}a_${gamma}g
    mkdir -p $res_dir
    cd $res_dir
    for randomness in 0.01 0.1 0.2 0.5 
    do
        for lambda in 0.0 0.1 0.5 0.9 0.99 #0.05 0.1 0.25 0.5 0.75 0.9 0.95
        do 
            echo "states" $n_states "lambda" $lambda "random" $randomness
	    time ${exec_dir}/online_algorithms $n_states $n_actions $gamma $lambda $randomness $n_runs $n_episodes $n_steps QLearning | grep REWARD >qlearning${lambda}l_${randomness}r.dat & 
	    time ${exec_dir}/online_algorithms $n_states $n_actions $gamma $lambda $randomness $n_runs $n_episodes $n_steps QLearningDirichlet | grep REWARD >qlearning_dir${lambda}l_${randomness}r.dat &
            wait
	done
    done
done