#! /bin/sh

mkdir -p $HOME/experiments
cd $HOME/experiments

script_dir=$HOME/projects/beliefbox/src/statistics/tests/scripts
environment=1DMaze

T=100000
iter=1000
maze_dir=$HOME/projects/beliefbox/dat
maze_name=maze1c
for action_rand in 1.0 0.5 0.1
do
    for maze_rand in 0 0.01 0.1
    do
        for n_obs in 2 16
        do
            out_dir=~/experiments/bvmdp/gridworld/${maze_name}/${n_obs}obs/p${maze_rand}_pa${action_rand}
            mkdir -p $out_dir
            echo -e "- action randomness:$action_rand\n- maze randomness: $maze_rand\n- observations: $n_obs\n- maze: $maze_name\n- iterations: $iter\n- T: $T" >${out_dir}/run.params
            
            for depth in 1 2 4 6 8
            do
                qsub -v"iter=$iter","T=$T","depth=$depth","maze_rand=$maze_rand","action_rand=$action_rand","maze_dir=$maze_dir","maze_name=$maze_name","n_obs=$n_obs","out_dir=$out_dir" ${script_dir}/factored_model_gridworld_sub.sh
            done
        done
    done
done

