#! /bin/bash

for gamma in 0.1 0.5 0.9 0.99 0.999 0.9999 0.99999 0.999999
do
    for method in 0 1
    do
        ./bin/toy_ucb_stopping_problem $method $gamma $iter 0 1000  >>quick_eval_ucb.out
        tail -n 1 quick_eval_ucb.out
    done
done
