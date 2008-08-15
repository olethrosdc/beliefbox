#! /bin/bash
#


results=results/belief_uct/bandit/uct_complete.out
experiments=1000
##echo "#meth #n_iter #gamma #actions #total #discounted #o_total #o_discounted" > $results

for actions in 2 4 8
do
    for iter in  1 2 3 4 5 6 7 8 10 12 14 18 22 28 32 64
    do
        for gamma in 0.5 0.9 0.99 0.999 
        do
            for method in 0 1 2 3 4 5 6 7 8 9 10 11 12 13
            do
                
                qsub -v"method=$method","iter=$iter","actions=$actions","experiments=$experiments"  ./bandit_uct_sub.sh
                tail -n 1 $results
            done
            
        done
    done
done