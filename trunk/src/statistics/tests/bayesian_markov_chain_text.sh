#! /bin/bash

horizon=
outdir=~/experiments/bpsr/alice
mkdir -p $outdir
summary=$outdir/summary.err
rm -f $summary
for (( i=0; i<=16; i++ ))
do
    ./bin/bayesian_markov_chain_text ./data/alice.dat $i $horizon | grep Err >>$summary
    mv bmc_text.error $outdir/bmc${i}.err
    mv bpsr_text.error $outdir/bpsr${i}.err
    mv ctw_text.error $outdir/ctw${i}.err
done

