#! /bin/sh


T=10000
iter=10
states=0
maze_dir=$HOME/projects/beliefbox/dat
maze_name=maze1c
for action_rand in 1.0 0.5 0.1
do
    for maze_rand in 0 0.01 0.1
    do
        for n_obs in 2 16
        do
            outdir=~/experiments/bvmdp/gridworld/${maze_name}/${n_obs}obs/p${maze_rand}_pa${action_rand}
            mkdir -p $outdir
            echo -e "- action randomness:$action_rand\n- maze randomness: $maze_rand\n- observations: $n_obs\n- maze: $maze_name\n- iterations: $iter\n- T: $T" >${outdir}/run.params
            
            for depth in 1 2 4 #6 8
            do
                for model in FMC BFMC BVMM CTW
                do
                    echo $model $depth
                    echo "./bin/factored_model Gridworld $model $iter $T $depth $maze_rand $action_rand $maze_dir/$maze_name $n_obs >${outdir}/${states}s_${depth}d_${model}.out"
                    time ./bin/factored_model Gridworld $model $iter $T $depth $maze_rand $action_rand $maze_dir/$maze_name $n_obs >${outdir}/${states}s_${depth}d_${model}.out
                done
            done
        done
    done
done
