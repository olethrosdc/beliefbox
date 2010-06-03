d=3;
n=10;
T=1000000;

for d in 3 2 4 5;
do
    outdir=$HOME/results/bvmdp/icml-workshop/${n}iter/Pendulum${d}
    mkdir -p $outdir
    
    nice -n 19 ./bin/pomdp_algorithms $d 3 0.99 0.8 0.01 $n 10000 $T QLearning Pendulum | grep PAYOFF >$outdir/QLearning.out &
    for i in 0 1 2
    do 
        nice -n 19 ./bin/pomdp_algorithms $d 3 0.99 0.8 0.01 $n 10000 $T BVMM Pendulum ${i} | grep PAYOFF >$outdir/vmdp${i}.out &
    done
done