#! /bin/bash

#data_name=dorian
#horizon=433520 # picture of dorian gray data size
#prior=0.02173913

#data_name=alice
#horizon=166862 # alice data size
#prior=0.02702703

data_name=splice
horizon=194590 # splice data size
#prior=0.1111111

prior=0.5

outdir=~/experiments/bpsr/${data_name}/p05
mkdir -p $outdir
cd $outdir
summary=$outdir/summary.err
summary_acc=$outdir/summary.acc
execdir=~/projects/beliefbox/src/statistics/tests/bin
datadir=~/projects/beliefbox/src/statistics/tests/data
rm -f $summary
rm -f $summary_acc
for (( i=0; i<=10; i++ ))
do
    ${execdir}/bayesian_markov_chain_text ${datadir}/${data_name}.txt $i $horizon $prior >out;
    grep Err out >>$summary
    grep Acc out >>$summary_acc
    tail $summary
    tail $summary_acc
    mv bmc_text.error bmc${i}.err
    mv bpsr_text.error bpsr${i}.err
    mv ctw_text.error ctw${i}.err
    mv polya_text.error polya${i}.err
    mv bmc_text.accuracy bmc${i}.acc
    mv bpsr_text.accuracy bpsr${i}.acc
    mv ctw_text.accuracy ctw${i}.acc
    mv polya_text.accuracy polya${i}.acc
done

