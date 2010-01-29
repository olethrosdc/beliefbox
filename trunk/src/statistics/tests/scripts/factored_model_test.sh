#! /bin/sh


iter=1000
T=100000
depth=0;

for ar in 1.0 0.5 0.1
do
    for r in 0.0 0.01 0.1
    do
        outdir=~/results/bvmdp/1DMaze/${r}p_${ar}ap
        mkdir -p $outdir
        echo "Params: $iter $T" >$outdir/params.txt
        for states in 4 6 8 10 12 16 32
        do
            for model in POMDP
            do
                echo "./bin/factored_model 1DMaze $model $iter $T $depth 0.0 1.0 $states >${outdir}/${states}s_${depth}d_${model}.out" >>$outdir/params.txt
                time ./bin/factored_model 1DMaze $model $iter $T 0 $r $ar $states  >${outdir}/${states}s_${depth}d_${model}.out
            done
#    for depth in 0 1 2 4 8 10 12
#    do
#        for model in FMC BFMC BVMM CTW
#        do
#            echo "./bin/factored_model 1DMaze $model $iter $T $depth 0.0 1.0 $s#tates >${outdir}/${states}s_${depth}d_${model}.out" >>$outdir/params.txt
#            time ./bin/factored_model 1DMaze $model $iter $T $depth 0.0 1.0 $states  >${outdir}/${states}s_${depth}d_${model}.out
#        done
        done
    done
done