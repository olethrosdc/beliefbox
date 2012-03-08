#!/bin/bash


N_STATES=10
N_RUNS=1
N_STEPS=10000
gamma=0.9
epsilon=0.5
alpha=0.1
lambda=0.0
samples=1
randomness=0

for algorithm in TdBma Sarsa QLearning # Model Sampling
do
    ./bin/online_algorithms --algorithm $algorithm --environment RandomMDP --maze_name $HOME/projects/beliefbox/dat/maze0 --n_runs $N_RUNS --n_steps $N_STEPS --alpha $alpha --epsilon $epsilon --gamma $gamma --lambda $lambda --max_samples $samples --randomness $randomness >$algorithm.out
    grep PAYOFF $algorithm.out > $algorithm.payoff
    grep EPISODE $algorithm.out > $algorithm.episode
    grep REWARD $algorithm.out > $algorithm.reward
done
