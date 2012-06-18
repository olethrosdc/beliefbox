runs=10000
steps=1000
gamma=0.95
epsilon=0.0

task=Chain

outdir=$HOME/results/$task/Fixed
mkdir -p $outdir


params="--environment $task --maze_name $HOME/projects/beliefbox/dat/$maze --gamma $gamma --epsilon 0.1 --n_runs $runs --n_steps $steps --n_episodes 10000 --episode_steps -1 --pit_value 0.0 --step_value 0.001"

for algorithm in Oracle QLearning Model UCRL
do
    echo $algorithm
    /usr/bin/time -o $outdir/${algorithm}.cpu ./bin/online_algorithms --algorithm $algorithm $params >$outdir/${algorithm}.out 
done

for algorithm in Oracle QLearning Model UCRL
do
    grep PAYOFF $outdir/${algorithm}.out > $outdir/${algorithm}.payoff
    grep RUN $outdir/${algorithm}.out > $outdir/${algorithm}.reward
    rm $outdir/${algorithm}.out
done


for prior in Fixed Normal Beta
do
    outdir=$HOME/results/$task/$prior
    mkdir -p $outdir
    
    params="--environment $task --maze_name $HOME/projects/beliefbox/dat/$maze --gamma $gamma --epsilon 0.0 --n_runs $runs --n_steps $steps --n_episodes 10000 --episode_steps -1 --pit_value 0.0 --step_value 0.001 --reward_prior $prior"



    for i in 001 002 004 008 016 032 064 128
    do 
        for algorithm in LSampling USampling
        do
            echo $algorithm $i; 
            /usr/bin/time -o $outdir/${algorithm}${i}.cpu  ./bin/online_algorithms --algorithm $algorithm --max_samples=${i} $params >$outdir/${algorithm}${i}.out &
        done 
        wait;
    done

    wait;

    for i in 001 002 004 008 016 032 064 128
    do 
        for algorithm in LSampling USampling
        do
            grep PAYOFF $outdir/${algorithm}${i}.out >$outdir/${algorithm}${i}.payoff; 
            grep RUN $outdir/${algorithm}${i}.out >$outdir/${algorithm}${i}.reward; 
            rm $outdir/${algorithm}${i}.out
        done
    done

done