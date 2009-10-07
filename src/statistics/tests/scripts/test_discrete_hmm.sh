#! /bin/bash

BINDIR=~/projects/beliefbox/src/statistics/tests/bin
RESDIR=~/experiments/discrete_hmm

mkdir -p $RESDIR

n_states=4
n_obs=4
stat=0.9
n_par=16
T=10000
n_iter=10

for n_obs in 2 4 8 16
do
    for n_states in  2 4 8 16
    do
        for n_par in 16
        do
            for stat in 0.5 0.9
            do
                resout=${n_states}s_${n_obs}o_${n_par}p_${stat}st_${T}T
                $BINDIR/discrete_hmm $n_states $n_obs $stat $n_par $T $n_iter >$RESDIR/${resout}.out &
            done
            wait;
            for stat in 0.5 0.9
            do
                resout=${n_states}s_${n_obs}o_${n_par}p_${stat}st_${T}T
                grep loss $RESDIR/${resout}.out   >$RESDIR/${resout}.loss
                grep pf_to_mix $RESDIR/${resout}.out   >$RESDIR/${resout}.cmp
                grep mixture $RESDIR/${resout}.out   >$RESDIR/${resout}.weights
            done
        done
    done
done
