experiments=100
rm -f quick_test.out
horizon=100000
for iter in 0 2 3 4 5
do
    for method in 0 1 2 3
    do
        time ./bin/bandit_uct_test $method 0.999999 $iter 2 0 $experiments $horizon >>quick_test.out
        tail quick_test.out
    done
done

### 5 iter 10^2 * 10^5 ~ 368 min -> 3.68 10^{-5} min per iteration

