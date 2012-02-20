runs=10000
steps=1000
gamma=0.95
epsilon=0.0

task=Chain
params="--n_states 5 --n_actions 2 --environment Chain --gamma $gamma --epsilon $epsilon --n_runs $runs --n_steps $steps"

outdir=$HOME/results/${task}/interval_100
mkdir -p $outdir
for algorithm in Model
do
    time ./bin/online_algorithms --algorithm $algorithm $params >$outdir/${algorithm}.out
    grep PAYOFF $outdir/${algorithm}.out > $outdir/${algorithm}.payoff
    grep RUN $outdir/${algorithm}.out > $outdir/${algorithm}.reward
    rm $outdir/${algorithm}.out
done

for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16;
do 
    echo $i; 
    time ./bin/online_algorithms --algorithm Sampling --max_samples=${i} $params >$outdir/sampling${i}.out &
done

wait

for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16;
do 
    grep PAYOFF $outdir/sampling${i}.out >$outdir/sampling${i}.payoff; 
    grep RUN $outdir/sampling${i}.out >$outdir/sampling${i}.reward; 
    rm $outdir/sampling${i}.out
done


