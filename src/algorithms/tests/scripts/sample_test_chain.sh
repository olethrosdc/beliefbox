runs=100
steps=1000
gamma=0.95
epsilon=0.0

task=Chain

outdir=$HOME/results/$task
mkdir -p $outdir


params="--environment $task --maze_name $HOME/projects/beliefbox/dat/$maze --gamma $gamma --epsilon 0.1 --n_runs $runs --n_steps $steps --n_episodes 10000 --episode_steps -1 --pit_value 0.0 --step_value 0.001"

for algorithm in Oracle QLearning Model UCRL
do
    /usr/bin/time -o $outdir/${algorithm}.cpu ./bin/online_algorithms --algorithm $algorithm $params >$outdir/${algorithm}.out 
    grep PAYOFF $outdir/${algorithm}.out > $outdir/${algorithm}.payoff
    grep RUN $outdir/${algorithm}.out > $outdir/${algorithm}.reward
    rm $outdir/${algorithm}.out
done

params="--environment $task --maze_name $HOME/projects/beliefbox/dat/$maze --gamma $gamma --epsilon 0.0 --n_runs $runs --n_steps $steps --n_episodes 10000 --episode_steps -1 --pit_value 0.0 --step_value 0.001"

for i in 1 2 4 6 8
do 
    echo $i; 
    /usr/bin/time -o $outdir/l_sampling${i}.cpu  ./bin/online_algorithms --algorithm Sampling --max_samples=${i} $params >$outdir/l_sampling${i}.out
done

wait;

for i in 1 2 4 6 8
do 
    grep PAYOFF $outdir/l_sampling${i}.out >$outdir/l_sampling${i}.payoff; 
    grep RUN $outdir/l_sampling${i}.out >$outdir/l_sampling${i}.reward; 
    rm $outdir/l_sampling${i}.out
done


for i in 1 2 4 6 8
do 
    echo $i; 
    /usr/bin/time -o $outdir/u_sampling${i}.cpu  ./bin/online_algorithms --algorithm Sampling --upper_bound --max_samples=${i} $params >$outdir/u_sampling${i}.out 
done

wait;

for i in 1 2 4 6 8
do 
    grep PAYOFF $outdir/u_sampling${i}.out >$outdir/u_sampling${i}.payoff; 
    grep RUN $outdir/u_sampling${i}.out >$outdir/u_sampling${i}.reward; 
    rm $outdir/u_sampling${i}.out
done

