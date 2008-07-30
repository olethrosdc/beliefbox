#! /bin/bash
#

for iter in  1 2 3 4 6 8 10 12 16 32 64 128
do
    results=results/bfs${iter}.out
    mkdir -p results
    qsub -v"iter=$iter","results=$results" ./test_toy_uct_sub.sh
done

