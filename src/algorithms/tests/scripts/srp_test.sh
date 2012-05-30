runs=1;
eps=1000;
maze=maze1
steps=100000
gamma=0.99;
for environment in Chain Optimistic ContextBandit;
do
    basedir=$HOME/results/srp

    for adversary in Random Fixed Heuristic;
    do
        resdir=$basedir/$environment/${adversary}
        mkdir -p $resdir
        for algo in Oracle LSampling USampling UCRL Model 
        do
            results=$resdir/$algo.out
            echo "$reward # $algo $adversary $environment";
            ./bin/srp --reward_prior Normal --maze_name ~/projects/beliefbox/dat/$maze --algorithm $algo --max_samples 4 --n_runs $runs --n_episodes $eps --environment $environment --adversary $adversary --epsilon 0.0 --gamma $gamma --n_states 16 >$results 

            grep REWARD $results >$resdir/$algo.reward
            grep EPISODE $results >$resdir/$algo.episode
            grep PAYOFF $results >$resdir/$algo.payoff
            rm $results
        done;
    done;
done
