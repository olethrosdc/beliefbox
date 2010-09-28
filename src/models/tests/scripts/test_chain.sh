#/ bin/bash

gamma=0.9
T=1000
samples=1000
for alpha in 0 0.5 1.0
do
    for knn in 3 5 10
    do
        for rbf in 0.1 1.0 10.0 100.0
        do
            echo "alpha" $alpha "knn" $knn "rbf" $rbf
            nice -n 19 ../bin/knn_test chain $alpha $knn $rbf $gamma $T $samples>out;
            grep V out >Vout;
            grep Q out>Qout;
            grep SAMP out >Sout;
            grep Adding out
            
            ##octave -q ./test_chain.m
            gnuplot test_chain.gnuplot
            outname=test_chain_${alpha}a_${knn}k_${rbf}r.eps
            mv test_chain.eps $outname
            #gv $outname&
        done
    done
done


