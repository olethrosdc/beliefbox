#! /bin/bash

samples=100
gamma=0.95
n_steps=10
accuracy=0.1
n_runs=1
task=Optimistic
randomness=1.0
epsilon=0.1

for n_demonstrations in 10
do
    for n_tasks in 2
    do 
        echo $n_tasks
        valgrind --leak-check=full ./bin/multitask_inverse_rl --MonteCarlo --n_chain_samples $samples --n_tasks $n_tasks --n_demonstrations $n_demonstrations --gamma $gamma --environment $task --randomness $randomness --n_steps $n_steps --accuracy $accuracy --n_runs $n_runs --epsilon $epsilon 2> valgrind.err
    done
done
