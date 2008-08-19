#! /bin/bash
#

for iter in  1 2 3 4 6 8 10 12 16 32 64 128
do
    results=results/bfs_slow_${iter}.out
    mkdir -p results
    rm -f $results
    bash ./test_toy_uct_sub2.sh $iter $results &
done

