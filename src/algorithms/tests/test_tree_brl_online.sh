#! /bin/bash

ulimit -Sv 1000000
n=100
for horizon in 0 1 2 3 4
do
	./bin/online_algorithms --environment Chain --n_states 5 --algorithm TBRL --horizon $horizon >out${horizon} 
done

