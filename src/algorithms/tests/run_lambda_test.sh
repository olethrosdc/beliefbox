#! /bin/bash

n_actions=2
gamma=0.9
exec_dir=~/projects/beliefbox/src/algorithms/tests/bin/

for n_states in 4 8 16 32 64 128
do

res_dir=~/experiments/prior/tests/${n_states}s_${n_actions}a_${gamma}g
mkdir -p $res_dir
cd $res_dir
for lambda in 0.0 0.5 0.9
do 
	for randomness in 0.0 0.1 0.5
	do
        echo "states" $n_states "lambda" $lambda "random" $randomness
	    ${exec_dir}/online_algorithms $n_states $n_actions $gamma $lambda $randomness | grep REWARD >out${lambda}l_${randomness}r.dat
	done
done
done