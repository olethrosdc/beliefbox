#! /bin/bash

samples=10000
gamma=0.95
n_steps=100
accuracy=0.001
n_runs=100
task=Optimistic
randomness=1.0
epsilon=0.1

for n_demonstrations in 100 200 500 1000
do
    for n_tasks in 1 2 5 10 20 50 100
    do 
        echo $n_tasks
        time ./bin/multitask_inverse_rl --MonteCarlo --n_chain_samples $samples --n_tasks $n_tasks --n_demonstrations $n_demonstrations --gamma $gamma --environment $task --randomness $randomness --n_steps $n_steps --accuracy $accuracy --n_runs $n_runs --epsilon $epsilon | tee opt_out_d${n_demonstrations}_t${n_tasks}
    done
done
