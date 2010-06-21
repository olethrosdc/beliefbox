#! /bin/bash

for gamma in 0.9999 #0.1 0.5 0.9 0.99 0.999
do
    for method in 0 1 2 3
    do
        for iter in 1 2 3 4 6 7 8 9 10 11 12 16 #1 2 3 4 5 6
        do

            ./bin/toy_uct_stopping_problem $method $gamma $iter 0 1000  >>quick_eval2.out
            tail quick_eval2.out
        done
    done
done