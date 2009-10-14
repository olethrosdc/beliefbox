#! /bin/bash

horizon=166862 # alice data size
for p in 0.5 0.001 0.01 0.1 0.9 0.99 0.999
  do 
  outdir=~/experiments/bpsr/alice/p${p}
  mkdir -p $outdir
  summary=$outdir/summary.err
  rm -f $summary
  for (( i=0; i<=16; i++ ))
    do
    ./bin/bayesian_markov_chain_text ./data/alice.txt $i $horizon $p| grep Err >>$summary
    tail $summary
    mv bmc_text.error $outdir/bmc${i}.err
    mv bpsr_text.error $outdir/bpsr${i}.err
    mv ctw_text.error $outdir/ctw${i}.err
    mv polya_text.error $outdir/polya${i}.err
  done
done

