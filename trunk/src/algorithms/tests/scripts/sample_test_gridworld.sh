runs=100
steps=10000
gamma=0.95
epsilon=0.1

#task=Chain
#params="--n_states 5 --n_actions 2 --environment Chain --gamma $gamma --epsilon $epsilon --n_runs $runs --n_steps $steps"

task=Gridworld


for randomness in 0.01 0.1 0.0
do
    for maze in maze0 maze1 maze2 # maze3 maze4
    do
        outdir=$HOME/results/$task/$maze/r${randomness}
        mkdir -p $outdir
        echo "writing resutls to $outdir"

        params="--environment $task --maze_name $HOME/projects/beliefbox/dat/$maze --gamma $gamma --n_runs $runs --n_steps $steps --n_episodes 10000 --episode_steps -1 --randomness $randomness --pit_value -1.0 --step_value 0.0"

        for algorithm in QLearning Model
        do
            /usr/bin/time -o $outdir/${algorithm}.cpu ./bin/online_algorithms --algorithm $algorithm $params --epsilon $epsilon >$outdir/${algorithm}.out &
        done
        wait;

        for algorithm in Oracle UCRL
        do
            /usr/bin/time -o $outdir/${algorithm}.cpu ./bin/online_algorithms --algorithm $algorithm $params --epsilon 0.0 >$outdir/${algorithm}.out &
        done
        wait;

        for algorithm in Oracle QLearning Model UCRL
        do
            grep PAYOFF $outdir/${algorithm}.out > $outdir/${algorithm}.payoff
            grep RUN $outdir/${algorithm}.out > $outdir/${algorithm}.reward
            rm $outdir/${algorithm}.out
        done

        params="--environment $task --maze_name $HOME/projects/beliefbox/dat/$maze --gamma $gamma --epsilon 0.0 --n_runs $runs --n_steps $steps --n_episodes 10000 --episode_steps -1 --randomness $randomness --pit_value -1.0 --step_value 0.01 --reward_prior Normal"


        for i in 001 002 004 008 016 032 064 128
        do 
            for algorithm in LSampling USampling
            do
                echo $algorithm $i; 
                /usr/bin/time -o $outdir/${algorithm}${i}.cpu  ./bin/online_algorithms --algorithm $algorithm --max_samples=${i} $params >$outdir/${algorithm}${i}.out &
            done 

            wait;

            for algorithm in LSampling USampling
            do
                grep PAYOFF $outdir/${algorithm}${i}.out >$outdir/${algorithm}${i}.payoff; 
                grep RUN $outdir/${algorithm}${i}.out >$outdir/${algorithm}${i}.reward; 
                rm $outdir/${algorithm}${i}.out
            done
        done


    done
done