
for gamma in 0.9 0.95 0.99 0.995 0.999
do
    for iter in 2 3 4 5 6 7 8 9 10
    do
        time ./bin/toy_uct_stopping_problem $gamma $iter 0 >>results
    done
done