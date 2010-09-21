model_name=BVMM
n_iter=10
T=10000
depth=1
for maze in maze1 maze6 #maze1 maze5
do
    for env_rand in 0.01 #0
    do
        for act_rand in 0.1 0.01 #0 0.1 #0.1 0.01 0
        do
            maze_dir=$HOME/projects/beliefbox/dat
            for observations in 2 16
            do 
                out_dir=$HOME/results/bvmdp/icml-workshop/${n_iter}iter/gridworld/${maze}/${observations}obs_${env_rand}er_${act_rand}ar/
                mkdir -p $out_dir
                echo "observations: " $observations "env rand: " $env_rand "act rand: " $act_rand "maze: " $maze >$out_dir/params
                for depth in 0 1 2 3 4
                do
                    echo "observations: " $observations "env rand: " $env_rand "act rand: " $act_rand "maze: " $maze "depth: " $depth
                    ./bin/factored_tree_rl Gridworld $model_name $n_iter $T $depth $env_rand $act_rand $maze_dir/$maze $observations >$out_dir/bvmdp_${depth}.out
    ##awk -v K=1000 -v col=2 -f ~/scripts/moving_average.awk out${depth} >out${depth}ma
                done
            done
        done
    done
done

##gdb --args ./bin/factored_tree_rl Gridworld $model_name $n_iter $T $depth $env_rand $act_rand $maze $observations
