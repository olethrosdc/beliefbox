runs=10
steps=10000
gamma=0.95
epsilon=0.0

## task=Chain
## params="--n_states 5 --n_actions 2 --environment Chain --gamma $gamma --epsilon $epsilon --n_runs $runs --n_steps $steps"

task=Gridworld
for maze in maze0 maze1 
do
    params="--environment $task --maze_name $HOME/projects/beliefbox/dat/$maze --gamma $gamma --epsilon $epsilon --n_runs $runs --n_steps $steps --n_episodes 10000 --episode_steps -1"
    
    outdir=$HOME/results/${task}/$maze
    mkdir -p $outdir
    for algorithm in QLearning Sarsa Model
    do
        time ./bin/online_algorithms --algorithm $algorithm $params >$outdir/${algorithm}.out
        grep PAYOFF $outdir/${algorithm}.out > $outdir/${algorithm}.payoff
        grep RUN $outdir/${algorithm}.out > $outdir/${algorithm}.reward
        rm $outdir/${algorithm}.out
    done
    
    for i in 1  2 3 4 5 6 7 8;
    do 
        echo $i; 
        time ./bin/online_algorithms --algorithm Sampling --max_samples=${i} $params >$outdir/sampling${i}.out;
        grep PAYOFF $outdir/sampling${i}.out >$outdir/sampling${i}.payoff; 
        grep RUN $outdir/sampling${i}.out >$outdir/sampling${i}.reward; 
        rm $outdir/sampling${i}.out
    done
done
