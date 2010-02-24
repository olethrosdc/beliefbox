#! /bin/bash

for dataset in cantrbry/alice29 cantrbry/asyoulik cantrbry/fields large/E.coli 
do
    datafile=~/data/${dataset}.dat
    outdir=~/experiments/bvmm/${dataset}
    mkdir -p $outdir
    for method in BMC BVMM CTW PPM;
    do
        echo "Method: $method Dataset $datafile"
        logfile=${outdir}/${method}.out;
        rm -f $logfile
        for depth in 1 2 3 4 5 6 7 8 10 12 14 16 20 24 28 32;
        do
            ./bin/context_tree_test $method $depth ${outdir}/${method}${depth}.out $datafile >>$logfile
            tail -n 1 $logfile
        done;
        grep loss $logfile | cut -f 3 -d ":" | cut -f 1 -d "," >${outdir}/${method}.loss
    done
done
