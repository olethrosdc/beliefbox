#! /bin/bash
#

results=results/belief_uct/bandit/baseline.out
experiments=1000

rm -f $results
for actions in 2 4 8 #16 32 64 128
do
    for gamma in  0.5 0.9 0.99 0.999 0.9999 0.99999 #0.99991 0.99992 0.99993 0.99994 0.99995 0.99996 0.99997 0.99998 0.99999 0.999991 0.999992 0.999993 0.999994 0.999995 0.999996 0.999997 0.999998 0.999999
    #for gamma in  0.5 0.9 0.99 0.999 0.9999 0.99999 0.999999
    do
        for method in 0 2 3 #4 5 6
        do
            time ./bin/bandit_ucb $method $gamma $actions 0 $experiments  >>$results 
            tail $results
        done
    done
done

# long only
#real    252m3.887s
#user    250m54.773s
#sys     0m25.562s

# long and double
# real    256m29.096s
# user    250m18.347s
# sys     0m25.890s
