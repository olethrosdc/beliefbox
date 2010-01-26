#! /bin/sh

outdir=~/experiments/bvmdp/1dmaze/2obs/
mkdir -p $outdir
iter=1000
T=100000
echo "Params: $iter $T" >$outdir/params.txt

for depth in 0 1 2 4 8 10 12 16
do
    for states in 1 2 3 4 5 6 7 8 16 32
    do
        for model in FMC BFMC BVMM CTW
        do
            time ./bin/factored_model 1DMaze $model $iter $T $depth 0.0 1.0 $states >${outdir}/${states}s_${depth}d_${model}.out
        done
    done
done