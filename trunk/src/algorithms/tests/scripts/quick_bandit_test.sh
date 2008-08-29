experiments=1000
horizon=100000
method=$1
for gamma in 0.5 0.9 0.99 0.999 0.9999 0.99999
do
    for iter in 0 1 2 3 4 5 8 #6 7 8 9 10 11 12
    do
        ./bin/bandit_uct $method $gamma $iter 2 0 $experiments >>uct_bnb.out
        tail uct_bnb.out
        done
    done
done

### 5 iter 10^2 * 10^5 ~ 368 min -> 3.68 10^{-5} min per iteration

