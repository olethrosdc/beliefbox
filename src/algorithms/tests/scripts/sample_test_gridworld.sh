runs=1
steps=10000
gamma=0.95
epsilon=0.0

## task=Chain
## params="--n_states 5 --n_actions 2 --environment Chain --gamma $gamma --epsilon $epsilon --n_runs $runs --n_steps $steps"

task=Gridworld


for randomness in 0.1
do
    for maze in maze1 maze2
    do
        outdir=$HOME/results/$task/$maze/r${randomness}
        mkdir -p $outdir

        params="--environment $task --maze_name $HOME/projects/beliefbox/dat/$maze --gamma $gamma --epsilon $epsilon --n_runs $runs --n_steps $steps --n_episodes 10000 --episode_steps -1 --randomness $randomness"

        for algorithm in Model
        do
            time ./bin/online_algorithms --algorithm $algorithm $params >$outdir/${algorithm}.out
            grep PAYOFF $outdir/${algorithm}.out > $outdir/${algorithm}.payoff
            grep RUN $outdir/${algorithm}.out > $outdir/${algorithm}.reward
            rm $outdir/${algorithm}.out
        done
        
        for i in 1 2 4 8;
        do 
            echo $i; 
            time ./bin/online_algorithms --algorithm Sampling --max_samples=${i} $params >$outdir/sampling${i}.out;
            grep PAYOFF $outdir/sampling${i}.out >$outdir/sampling${i}.payoff; 
            grep RUN $outdir/sampling${i}.out >$outdir/sampling${i}.reward; 
            rm $outdir/sampling${i}.out
        done
    done
done