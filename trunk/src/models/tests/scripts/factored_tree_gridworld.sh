model_name=BVMM
n_iter=10
T=100000
depth=1
env_rand=0.01
act_rand=0.1
maze=~/projects/beliefbox/dat/maze1
observations=2

for depth in 1 2 4 8 16; 
do
    ./bin/factored_tree_rl Gridworld $model_name $n_iter $T $depth $env_rand $act_rand $maze $observations >out${depth}
    awk -v K=1000 -v col=2 -f ~/scripts/moving_average.awk out${depth} >out${depth}ma
done


##gdb --args ./bin/factored_tree_rl Gridworld $model_name $n_iter $T $depth $env_rand $act_rand $maze $observations
