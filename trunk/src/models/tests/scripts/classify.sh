datadir=$HOME/data/classification
resdir=$HOME/results/classification
mkdir -p $resdir
for dataset in flare iris optdigits pendigits spect
do
    train=$datadir/$dataset/train.dat
    test=$datadir/$dataset/test.dat
    echo ""
    echo "Data: $dataset"
    echo "Training set: $train"
    echo "Test set: $test"
    outdir=$resdir/$dataset
    mkdir -p $outdir
	for knn in 1 2 3
	do
#		fname=$outdir/knn${knn}
#		/usr/bin/time --verbose --output $fname.time ./bin/classify --train $train --test $test --randomise --knn $knn  >$fname.out
#		grep TRAIN $fname.out >${fname}_train.res
#		grep TEST $fname.out >${fname}_test.res
        
	    for tree in 32 64 128 #1 2 4 8 16
	    do
		    fname=$outdir/tree${tree}_knn${knn}
		    /usr/bin/time --verbose --output $fname.time ./bin/classify --train $train --test $test --randomise --tree $tree --knn $knn >$fname.out
		    grep TRAIN $fname.out >${fname}_train.res
		    grep TEST $fname.out >${fname}_test.res
        done
    done

#		fname=$outdir/gauss
#		/usr/bin/time --verbose --output $fname.time ./bin/classify --train $train --test $test --randomise --gaussian >$fname.out
#		grep TRAIN $fname.out >${fname}_train.res
#		grep TEST $fname.out >${fname}_test.res
        
	    for tree in 32 64 128 #1 2 4 8 16
	    do
		    fname=$outdir/tree${tree}_gauss
		    /usr/bin/time --verbose --output $fname.time ./bin/classify --train $train --test $test --randomise --tree $tree --gaussian >$fname.out
		    grep TRAIN $fname.out >${fname}_train.res
		    grep TEST $fname.out >${fname}_test.res
        done


    echo "---------------------------------------------------------------------"
done