# shell for the job:
#PBS -S /bin/bash
# use one node with 8 cores:
#PBS -lnodes=1:ppn=8
# job requires at most 16 hours, 0 minutes
#     and 0 seconds wallclock time
#PBS -lwalltime=16:00:00
# cd to the directory where the program is to be called:

runs=10000
steps=1000
gamma=0.95
epsilon=0.0


task=Optimistic
params="--n_states 5 --n_actions 2 --environment Optimistic --gamma $gamma --epsilon $epsilon --n_runs $runs --n_steps $steps"

## task=Chain
## params="--n_states 5 --n_actions 2 --environment Chain --gamma $gamma --epsilon $epsilon --n_runs $runs --n_steps $steps"

## task=RiverSwim
## params="--environment RiverSwim --gamma $gamma --epsilon $epsilon --n_runs $runs --n_steps $steps"

outdir=$HOME/results/${task}/epsilon${epsilon}
mkdir -p $outdir
#for algorithm in QLearning Sarsa Model
#do
#    time ./bin/online_algorithms --algorithm $algorithm $params >$outdir/${algorithm}.out
#    grep PAYOFF $outdir/${algorithm}.out > $outdir/${algorithm}.payoff
#    grep RUN $outdir/${algorithm}.out > $outdir/${algorithm}.reward
#    rm $outdir/${algorithm}.out
#done

for i in 1 2 3 4 5 6 7 8;
do 
    echo $i; 
    time ./bin/online_algorithms --algorithm Sampling --max_samples=${i} $params >$outdir/sampling${i}.out; 
    grep PAYOFF $outdir/sampling${i}.out >$outdir/sampling${i}.payoff; 
    grep RUN $outdir/sampling${i}.out >$outdir/sampling${i}.reward; 
    rm $outdir/sampling${i}.out
done
