#! /bin/bash

n_iter=100
n_states=8
n_observations=8
stationarity=0.9
T=100
out_dir=~/experiments/hmm_em/${n_states}s_${n_observations}s_${stationarity}st_${T}T_${n_iter}iter
mkdir -p $out_dir
for n_em_iter in 1 2 3 4 5 6 7 8 9 10 12 14 16 20 25
do
    ./bin/test_hmm_em $n_states $n_observations $stationarity $n_em_iter $T $n_iter >$out_dir/${n_em_iter}em.out
done
