#! /bin/bash

n_actions=4
gamma=0.9
exec_dir=~/projects/beliefbox/src/algorithms/tests/bin/

for n_states in 4 16 64 256
do
    res_dir=~/experiments/prior/tests/${n_states}s_${n_actions}a/
    mkdir -p $res_dir
    cd $res_dir
    for lambda in 0.0 0.5 0.9
    do 
	for randomness in 0.0 0.1 0.5
	do
	    ${exec_dir}/online_algorithms $n_states $n_actions $gamma $lambda $randomness | grep REWARD >out${lambda}l_${randomness}r.dat
	done
    done
done