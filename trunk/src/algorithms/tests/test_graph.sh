#! /bin/bash

gamma=0.99
depth=4
horizon=10
./bin/toy_uct_stopping_problem $gamma $depth 0 1 $horizon

for (( i=0; i<horizon; i++ ));
do
    dot -Tps -o test${i}.ps test${i}.dot
done
