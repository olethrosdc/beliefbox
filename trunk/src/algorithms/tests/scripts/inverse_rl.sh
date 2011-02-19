baseres=$HOME/results/inverse_rl
n_runs=1000

for iter in 1 10 100 1000 10000
do

    environment=Optimistic
    resdir=$baseres/$environment
    mkdir -p $resdir
    echo $environment $n_states
    ./bin/inverse_rl  --environment $environment --n_runs $n_runs --iterations $iter | grep DV >$resdir/${iter}i_DV.out

    environment=Gridworld
    for maze in maze0 maze1 maze2 maze3
    do
        resdir=$baseres/$environment/$maze
        mkdir -p $resdir
        for n_states in 2 4 8 16 32 64 128
        do
            echo $environment $maze $n_states
            ./bin/inverse_rl  --environment $environment --maze_name $maze --n_runs $n_runs --iterations $iter | grep DV >$resdir/${iter}i_DV.out
        done
    done

    environment=Chain
    resdir=$baseres/$environment
    mkdir -p $resdir
    for n_states in 2 4 8 16 #32 64 128
    do
        echo $environment $n_states
        ./bin/inverse_rl --n_states $n_states --environment $environment --n_runs $n_runs --iterations $iter | grep DV >$resdir/${n_states}s_${iter}i_DV.out
    done

    environment=RandomMDP
    resdir=$baseres/$environment
    mkdir -p $resdir
    for n_states in 4 8 16 32 64 128 256
    do
        echo $environment $n_states
        ./bin/inverse_rl --n_states $n_states --environment $environment --n_runs $n_runs --iterations $iter | grep DV >$resdir/${n_states}s_${iter}i_DV.out
    done

done

