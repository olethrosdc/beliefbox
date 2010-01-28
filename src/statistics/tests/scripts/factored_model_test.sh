#! /bin/sh

outdir=~/experiments/bvmdp/1dmaze/2obs/
mkdir -p $outdir
iter=100
T=100000
echo "Params: $iter $T" >$outdir/params.txt

for states in 4 6 8 10 12 16 32
do
    for depth in 0 1 2 4 8 10 12
    do
        for model in FMC BFMC BVMM CTW
        do
            echo "./bin/factored_model 1DMaze $model $iter $T $depth 0.0 1.0 $states >${outdir}/${states}s_${depth}d_${model}.out" >>$outdir/params.txt
            time ./bin/factored_model 1DMaze $model $iter $T $depth 0.0 1.0 $states  >${outdir}/${states}s_${depth}d_${model}.out
        done
    done
done