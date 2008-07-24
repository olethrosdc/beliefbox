#! /bin/bash

rm -f results
for gamma in 0.9 0.91 0.925 0.95 0.99 0.995 0.999
do
    for iter in 2 4 8 16 32 64
    do
        time ./bin/toy_uct_stopping_problem $gamma $iter 0 >>results
        tail -n 1 results
    done
done