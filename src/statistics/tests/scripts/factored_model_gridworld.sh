#! /bin/sh

outdir=~/experiments/bvmdp/gridworld/maze1
mkdir -p $outdir
T=100000
iter=1000
states=0
echo "factored_model $iter $T" >${outdir}/run.params

for depth in 0 1 2 3 4 5 6 7 8 9 10
do
    for model in FMC BFMC BVMM CTW
    do
        time ./bin/factored_model $iter $T $states $depth $model >${outdir}/${states}s_${depth}d_${model}.out
    done
done