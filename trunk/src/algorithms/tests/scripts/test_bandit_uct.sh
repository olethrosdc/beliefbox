#! /bin/bash
#


results=results/belief_uct/bandit/uct_is_q.out
experiments=100
#echo "#meth #n_iter #gamma #actions #total #discounted #o_total #o_discounted" > $results

#rm -f $results
for actions in 2 4 #8
do
    for gamma in  0.5 0.9 0.99 0.999 0.9999 0.99999
    do
        
        for iter in  1 2 3 4 5 6 7 8
        do
            for method in 0 1 4 5 13 15 16 17
	        do
	        #qsub -v"method=$method","iter=$iter","actions=$actions","experiments=$experiments","results=$results"  ./test_bandit_uct_sub.sh
                time ./bin/bandit_uct $method $gamma $iter $actions 0 $experiments  >>$results 
            done
        done
    done
done