#! /bin/bash

for gamma in 0.1 0.5 0.9 0.99 0.999 #0.1 0.5 0.9 0.99 0.999
do
    for method in 4 5 6 7
    do
        for iter in 1 2 3 4 5 6 7 8 9 10 11 12 16
        do
            ./bin/toy_uct_stopping_problem $method $gamma $iter 0 1000  >>quick_eval3.out
            tail -n 1 quick_eval3.out
        done
    done
done