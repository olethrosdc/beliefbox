#! /bin/sh

for p in 0.0 0.01 0.1
do
    outdir=~/experiments/bvmdp/1dmaze/1000iter/2obs/p${p}
    mkdir -p $outdir
    iter=1000
    T=100000
    echo "Params: $iter $T" >$outdir/params.txt
    
    for states in 4 8 10 12 16 32
    do
        for depth in 1 2 4 8 10 12
        do
            for model in FMC BFMC BVMM CTW
            do
                echo "./bin/factored_model 1DMaze $model $iter $T $depth $p 1.0 $states >${outdir}/${states}s_${depth}d_${model}.out" >>$outdir/params.txt
                time ./bin/factored_model 1DMaze $model $iter $T $depth $p 1.0 $states  >${outdir}/${states}s_${depth}d_${model}.out
            done
        done
    done
done