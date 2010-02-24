#! /bin/bash

for method in BMC BVMM CTW PPM;
do
    rm -f ${method}.out
    for depth in 1 2 3 4 5 6 7 8;
    do
        ./bin/context_tree_test $method $depth ${method}${depth}.out data/dorian.txt >>${method}.out;
        tail -n 1 ${method}.out
    done;
    grep loss ${method}.out | cut -f 3 -d ":" | cut -f 1 -d "," >${method}.loss
done
