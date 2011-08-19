#! /bin/bash

samples=1000
gamma=0.95
n_steps=100
accuracy=0.1
n_runs=10
task=Optimistic
randomness=1.0
n_demonstrations=100
epsilon=0.1

for n_tasks in 1 2 5 10 20 50 100
do 
    echo $n_tasks
    time ./bin/multitask_inverse_rl --MonteCarlo --n_chain_samples $samples --n_tasks $n_tasks --n_demonstrations $n_demonstrations --gamma $gamma --environment $task --randomness $randomness --n_steps $n_steps --accuracy $accuracy --n_runs $n_runs --epsilon $epsilon | tee out${n_tasks}
done
