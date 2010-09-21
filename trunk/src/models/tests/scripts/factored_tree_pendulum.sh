model_name=BVMM
n_iter=10
T=100000
depth=1
env_rand=0.01;
environment=Pendulum
for act_rand in 0.0 #0 0.1 #0.1 0.01 0
do
    maze_dir=$HOME/projects/beliefbox/dat
    for observations in 4 2 6 5 3
    do 
        out_dir=$HOME/results/bvmdp/icml-workshop/${n_iter}iter/${environment}/g0.99_${observations}obs_${env_rand}er_${act_rand}ar/
        mkdir -p $out_dir
        echo "observations: " $observations "env rand: " $env_rand "act rand: " $act_rand "T: " $T >$out_dir/params
        for depth in 0;# 1 2 3 4;
        do
            echo "observations: " $observations "env rand: " $env_rand "act rand: " $act_rand "T: " $T "depth: " $depth
            nice -n 19 ./bin/factored_tree_rl ${environment} $model_name $n_iter $T $depth $env_rand $act_rand $observations  >$out_dir/bvmdp_${depth}.out 
    ##awk -v K=1000 -v col=2 -f ~/scripts/moving_average.awk out${depth} >out${depth}ma
        done
        wait;
    done
done

##gdb --args ./bin/factored_tree_rl Gridworld $model_name $n_iter $T $depth $env_rand $act_rand $maze $observations
