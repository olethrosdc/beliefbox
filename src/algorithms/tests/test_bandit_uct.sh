#! /bin/bash
#

actions=2
results=results/belief_uct/bandit/uct.out
experiments=1000
echo "#meth #n_iter #gamma #actions #total #discounted #o_total #o_discounted" > $results

for iter in  1 2 3 4 5 6 7 8 10 12 14 18 22
do
    for method in 0 1 2 3 4
    do
        for gamma in 0.1 0.5 0.9 0.99
        do
            ./bin/bandit_uct $method $gamma $iter $actions 0 $experiments >>$results
            tail -n 1 $results
        done
        
    done
done

