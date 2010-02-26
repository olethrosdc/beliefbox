#! /bin/bash

#for dataset in calgary/trans calgary/paper1 calgary/pic #calgary/bib calgary/book1 calgary/book2 calgary/geo
for dataset in large/E.coli calgary/bib calgary/book1 calgary/book2 calgary/geo calgary/news calgary/obj1 calgary/obj2 calgary/paper1 calgary/paper2 calgary/pic calgary/progc calgary/progl calgary/progp calgary/trans
do
    datafile=~/data/${dataset}.dat
    outdir=~/experiments/bvmm/${dataset}
    mkdir -p $outdir
    for method in  PPM BMC BVMM #CTW PPM;
    do
        echo "Method: $method Dataset $datafile"
        logfile=${outdir}/${method}.out;
        rm -f $logfile
        for depth in 1 2 3 4 5 6 7 8 9 10 11 12  14 16 20 24 28 32;
        do
            ./bin/context_tree_test $method $depth ${outdir}/${method}${depth}.out $datafile >>$logfile
            tail -n 1 $logfile
        done;
        grep loss $logfile | cut -f 3 -d ":" | cut -f 1 -d "," >${outdir}/${method}.loss
    done
done
