#! /bin/bash
#



#echo "#meth #n_iter #gamma #actions #total #discounted #o_total #o_discounted" > $results
experiments=1000

for gamma in  0.5 #0.9 0.99 0.999 0.9999 0.99999
  do
  for iter in  0 1 #2 3 4 5 6 7 8 10 12 14 18 22 28 32 64
    do
    for method in 0 1 2 3 4 5 6 7 8 9 10 11 12 13
      do
      results=results/belief_uct/bandit/g${gamma}_m${method}_i${iter}_uct_complete.out
      qsub -v"method=$method","iter=$iter","gamma=$gamma","experiments=$experiments","results=$results"  ./test_bandit_uct_sub.sh
    done
  done
done