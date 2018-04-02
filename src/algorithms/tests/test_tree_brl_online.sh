#! /bin/bash

ulimit -Sv 2000000
n=100
for horizon in 0 1 2 3 4
do
	echo $horizon
	./bin/online_algorithms --environment Chain --n_states 5 --n_episodes 1 --n_runs 1 --gamma 0.9 --n_steps 10000 --algorithm TBRL --horizon $horizon 
done

