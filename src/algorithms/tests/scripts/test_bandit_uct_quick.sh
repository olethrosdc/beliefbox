#! /bin/bash
#


results=results/belief_uct/bandit/uct_quick.out
experiments=1000
##echo "#meth #n_iter #gamma #actions #total #discounted #o_total #o_discounted" > $results

for actions in 2 4 #8
do
    for iter in  0 2 5 21
    do
        for method in 0 1 3 9
	    do
	        #qsub -v"method=$method","iter=$iter","actions=$actions","experiments=$experiments","results=$results"  ./test_bandit_uct_sub.sh
            for gamma in  0.5 0.9 0.99 0.999 0.9999 0.99991 0.99992 0.99993 0.99994 0.99995 0.99996 0.99997 0.99998 0.99999 0.999991 0.999992 0.999993 0.999994 0.999995 0.999996 0.999997 0.999998 0.999999
            do
                time ./bin/bandit_uct $method $gamma $iter $actions 0 $experiments  >>$results 
            done
        done
    done
done