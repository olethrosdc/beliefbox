for a in 2 4 8;
do
    for gamma in 0.99 0.999
    do
        for i in 0 1 2 4 8
        do
            for m in 0 1 2 3 4 5 6 13 14 15
            do
                nice ./bin/bandit_uct_rebuild_tree $m $gamma $i $a 0 100 >>benchmark_rebuild_tree.out;
                nice ./bin/bandit_uct $m $gamma $i $a 0 100 >>benchmark_keep_treeout;
            done
        done
    done
done
