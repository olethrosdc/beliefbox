# shell for the job:
#PBS -S /bin/bash
# use one node with 16 cores:
#PBS -lnodes=1:ppn=16
# job requires at most 1 hours, 0 minutes
#     and 0 seconds wallclock time
#PBS -lwalltime=1:00:00
# cd to the directory where the program is to be called:

resdir=$basedir/$environment
mkdir -p $resdir
do
    for adversary in Random Fixed Heuristic;
    do
        for algo in Oracle LSampling USampling UCRL Model 
        do
            results=$resdir/$algo.out
            echo "$reward # $algo $adversary $environment";
            ./bin/srp --reward_prior Normal --maze_name ~/projects/beliefbox/dat/$maze --algorithm $algo --max_samples 4 --n_runs $runs --n_episodes $eps --environment $environment --adversary $adversary --epsilon 0.0 >$results &
        done;
    done;
done
wait;

