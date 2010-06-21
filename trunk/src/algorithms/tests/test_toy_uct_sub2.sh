#! /bin/bash
#
iter=$1
results=$2

for gamma in  0.1 0.5 0.9 0.99 0.999 0.9999 0.99999 0.999999
do
    time nice -n 19 ./bin/toy_uct_stopping_problem $gamma $iter 0 10000 >>$results 
done

