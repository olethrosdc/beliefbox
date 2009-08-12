#! /bin/bash
outdir=~/experiments/bpsr/splice/
mkdir -p $outdir
rm -f $outdir/summary.err
horizon=100000
for (( i=0; i<=16; i++ ))
do
    ./bin/bayesian_markov_chain_text ./data/splice.txt $i $horizon | grep Err  >>$outdir/summary.err
    mv bmc_text.error $outdir/bmc${i}.err
    mv bpsr_text.error $outdir/bpsr${i}.err
    mv ctw_text.error $outdir/ctw${i}.err
done

