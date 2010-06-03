d=3;
n=1000;
T=1000000;

for d in 3 2 4 5;
do
    outdir=$HOME/results/bvmdp/icml-workshop/${n}iter/Pendulum${d}
    mkdir -p $outdir
    
    for i in 3 4 8 16
    do 
        nice -n 19 ./bin/pomdp_algorithms $d 3 0.99 0.8 0.01 $n 10000 $T BVMM Pendulum ${i} | grep PAYOFF >$outdir/vmdp${i}.out &
    done
    wait
done