#! /bin/bash
#


results=results/belief_uct/bandit/uct2.out
experiments=1000
##echo "#meth #n_iter #gamma #actions #total #discounted #o_total #o_discounted" > $results

for actions in 4 8
do
    for iter in  2 3 4 5 6 7 8 10 12 14 18 22
    do
        for gamma in 0.1 0.5 0.9 0.99 #0.999
        do
            for method in 12 #8 9 10 11 # 0 1 2 3 4 5 6 7
            do
                ./bin/bandit_uct $method $gamma $iter $actions 0 $experiments  >>$results
                tail -n 1 $results
            done
            
        done
    done
done