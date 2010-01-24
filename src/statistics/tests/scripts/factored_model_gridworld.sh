#! /bin/sh

outdir=~/experiments/bvmdp/gridworld/maze1
mkdir -p $outdir
T=10000
states=0
for depth in 0 1 2 4 6 8 10
do
    for model in FMC BVMM CTW
    do
        time ./bin/factored_model $T $states $depth $model >${outdir}/${states}s_${depth}d_${model}.out
    done
done