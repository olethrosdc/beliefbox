#! /bin/sh

outdir=~/experiments/bvmdp/1dmaze
mkdir -p $outdir
T=100000
for depth in 0 1 2 4 8 16
do
    for states in 1 2 3 4 5 6 7 8 16 32
    do
        for model in FMC BVMM CTW
        do
            time ./bin/factored_model $T $states $depth $model >${outdir}/${states}s_${depth}d_${model}.out
        done
    done
done