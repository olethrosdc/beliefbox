#! /bin/bash

for n_actions in 2 4 8 16 32 64 128
do
    for gamma in 0.1 0.5 0.9 0.99 0.999 0.9999 0.99999 0.999999
    do
        for method in 0 1
        do
            ./bin/bandit_ucb $method $gamma $n_actions 0 1000  >>bandit_eval_ucb.out
            tail -n 1 bandit_eval_ucb.out
        done
    done
done

