#! /bin/bash

data_name=dorian
horizon=433520 # picture of dorian gray data size
prior=0.02173913

#data_name=alice
#horizon=166862 # alice data size
#prior=0.02702703

#data_name=splice
#horizon=194590 # splice data size
#0.1111111

outdir=~/experiments/bpsr/${data_name}/KpowN
mkdir -p $outdir
summary=$outdir/summary.err
rm -f $summary
for (( i=0; i<=16; i++ ))
do
    ./bin/bayesian_markov_chain_text ./data/${data_name}.txt $i $horizon $prior | grep Err >>$summary
    tail $summary
    mv bmc_text.error $outdir/bmc${i}.err
    mv bpsr_text.error $outdir/bpsr${i}.err
    mv ctw_text.error $outdir/ctw${i}.err
    mv polya_text.error $outdir/polya${i}.err
    mv bmc_text.accuracy $outdir/bmc${i}.acc
    mv bpsr_text.accuracy $outdir/bpsr${i}.acc
    mv ctw_text.accuracy $outdir/ctw${i}.acc
    mv polya_text.accuracy $outdir/polya${i}.acc
done

