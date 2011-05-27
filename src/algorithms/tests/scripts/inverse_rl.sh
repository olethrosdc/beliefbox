baseres=$HOME/results/inverse_rl
n_runs=100

for steps in 100 200 400 800 1600 3200 6400 12800 25600 51200 102400
do

    environment=Optimistic
    resdir=$baseres/$environment
    mkdir -p $resdir
    for gamma in 0.9 0.95 0.99 0.999
    do
        echo $environment $gamma
        ./bin/inverse_rl  --environment $environment --gamma $gamma --n_runs $n_runs --n_steps $steps | tee $resdir/${steps}T_${gamma}g.out | grep DV >$resdir/${steps}T_${gamma}g_DV.out &
    done
    wait;

    environment=Gridworld
    for maze in maze0 maze1 maze2 maze3
    do
        resdir=$baseres/$environment/$maze
        mkdir -p $resdir
        echo $environment $maze
        maze_path=$HOME/projects/beliefbox/dat/$maze
        ./bin/inverse_rl  --environment $environment --maze_name $maze_path --n_runs $n_runs --n_steps $steps | tee $resdir/${steps}T.out | grep DV >$resdir/${steps}T_DV.out &
    done
    wait;

    environment=Chain
    resdir=$baseres/$environment
    mkdir -p $resdir
    for n_states in 2 4 8 16 
    do
        echo $environment $n_states
        ./bin/inverse_rl --n_states $n_states --environment $environment --n_runs $n_runs --n_steps $steps | tee $resdir/${n_states}s_${steps}T.out | grep DV >$resdir/${n_states}s_${steps}T_DV.out &
    done
    wait;

    environment=RiverSwim
    resdir=$baseres/$environment
    mkdir -p $resdir
    for n_states in 3 4 5 6
    do
        echo $environment $n_states
        ./bin/inverse_rl --n_states $n_states --environment $environment --n_runs $n_runs --n_steps $steps | tee $resdir/${n_states}s_${steps}T.out | grep DV >$resdir/${n_states}s_${steps}T_DV.out &
    done
    wait;

    environment=RandomMDP
    resdir=$baseres/$environment
    mkdir -p $resdir
    for n_states in 4 16 64 256
    do
        echo $environment $n_states
        ./bin/inverse_rl --n_states $n_states --environment $environment --n_runs $n_runs --n_steps $steps | tee >$resdir/${n_states}s_${steps}T.out | grep DV >$resdir/${n_states}s_${steps}T_DV.out &
    done
    wait;

done

