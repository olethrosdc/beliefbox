#! /bin/bash

samples=100;#00
gamma=0.95
n_steps=100
accuracy=0.01
n_runs=2;#100
task=Gridworld
randomness=0.1
epsilon=0.1
maze=$HOME/projects/beliefbox/dat/maze2
for n_demonstrations in 10 #20 50 100
do
    for n_tasks in 2;#1 2 5 10
    do 
        echo $n_tasks
        valgrind ./bin/multitask_inverse_rl --MonteCarlo --n_chain_samples $samples --n_tasks $n_tasks --n_demonstrations $n_demonstrations --gamma $gamma --environment $task --randomness $randomness --n_steps $n_steps --accuracy $accuracy --n_runs $n_runs --epsilon $epsilon --maze_name $maze 2>valgrind.err;# | tee maze_out_d${n_demonstrations}_t${n_tasks}
    done
done
