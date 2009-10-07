#! /bin/bash

exec_dir=$HOME/projects/beliefbox/src/statistics/tests/bin

iter=100
max_states=4
for T in 100 1000 10000
  do
  for n_states in 2 4 8
    do
    for n_obs in 2 4 8
      do
      for max_states in 4 8
        do
        out_dir=~/experiments/bpsr/T${T}_ms${max_states}_o${n_obs}_s${n_states}_iter${iter}
        mkdir -p $out_dir
        cd $out_dir
        ${exec_dir}/bayesian_markov_chain $T $iter $n_obs $max_states $n_states &
      done
      wait;
    done
  done
done
