baseres=$HOME/results/inverse_rl
n_runs=1000

for iter in 1 10 100 1000 10000
do

    environment=Optimistic
    resdir=$baseres/$environment
    mkdir -p $resdir
    for gamma in 0.9 0.95 0.99 0.999
    do
        echo $environment $gamma
        ./bin/inverse_rl  --environment $environment --gamma $gamma --n_runs $n_runs --iterations $iter | grep DV >$resdir/${iter}i_${gamma}g_DV.out &
    done
    wait;

    environment=Gridworld
    for maze in maze0 maze1 maze2 maze3
    do
        resdir=$baseres/$environment/$maze
        mkdir -p $resdir
        echo $environment $maze
        maze_path=$HOME/projects/beliefbox/dat/$maze
        ./bin/inverse_rl  --environment $environment --maze_name $maze_path --n_runs $n_runs --iterations $iter | grep DV >$resdir/${iter}i_DV.out &
    done
    wait;

    environment=Chain
    resdir=$baseres/$environment
    mkdir -p $resdir
    for n_states in 2 4 8 16 
    do
        echo $environment $n_states
        ./bin/inverse_rl --n_states $n_states --environment $environment --n_runs $n_runs --iterations $iter | grep DV >$resdir/${n_states}s_${iter}i_DV.out &
    done
    wait;

    environment=RandomMDP
    resdir=$baseres/$environment
    mkdir -p $resdir
    for n_states in 4 16 64 256
    do
        echo $environment $n_states
        ./bin/inverse_rl --n_states $n_states --environment $environment --n_runs $n_runs --iterations $iter | grep DV >$resdir/${n_states}s_${iter}i_DV.out &
    done
    wait;

done

