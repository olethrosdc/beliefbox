#! /bin/bash

topdir=~/experiments/bpsr/hmm_mc

T=10000
n_iter=100
n_obs=4
n_hist=4
n_mc_states=4

for n_obs in 2 4 8
do
    for n_mc_states in 2 4 8
    do
        for n_hist in 1 2 3 4 5 6 7 8 10 12 16
        do
            resdir=$topdir/${T}T_${n_iter}iter_${n_obs}o_${n_hist}h_${n_mc_states}s
            mkdir -p $resdir
            cd $resdir
            pwd
            ~/projects/beliefbox/src/statistics/tests/bin/bayesian_markov_chain $T $n_iter $n_obs $n_hist $n_mc_states
        done
    done
done

