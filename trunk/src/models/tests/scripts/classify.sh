datadir=~/data/classification
resdir=~/results/classification
mkdir -p $resdir
for dataset in iris optdigits pendigits spect
do
    train=$datadir/$dataset/train.dat
    test=$datadir/$dataset/test.dat
    echo ""
    echo "Data: $dataset"
    echo "Training set: $train"
    echo "Test set: $test"
    outdir=$resdir/$dataset
    mkdir -p $outdir
#	for knn in 1 2 3 4 5 6
#	do
#		fname=$outdir/knn${knn}
#		/usr/bin/time --verbose --output $fname.time ./bin/classify --train $train --test $test --randomise --knn $knn  >$fname.out
#		grep TRAIN $fname.out >${fname}_train.res
#		grep TEST $fname.out >${fname}_test.res
#		for mix in 1 2 4 8 16
#		do
#			fname=$outdir/knn${knn}_mix${mix}
#			/usr/bin/time --verbose --output $fname.time ./bin/classify --train $train --test $test --randomise --knn $knn --mixture $mix >$fname.out
#			grep TRAIN $fname.out >${fname}_train.res
#			grep TEST $fname.out >${fname}_test.res
#		done
#	done

	fname=$outdir/tree${knn}
	for tree in 1 2 4 8 16 32 64 128
	do
		fname=$outdir/tree${tree}
		/usr/bin/time --verbose --output $fname.time ./bin/classify --train $train --test $test --randomise --tree $tree >$fname.out
		grep TRAIN $fname.out >${fname}_train.res
		grep TEST $fname.out >${fname}_test.res
    done

# fname=$outdir/sparse${ktreenn}
#	for sparse in 1 2 4 8 16 32 64 128
#	do
#		fname=$outdir/sparse${sparse}
#		/usr/bin/time --verbose --output $fname.time ./bin/classify --train $train --test $test --randomise --sparse $sparse >$fname.out
#		grep TRAIN $fname.out >${fname}_train.res
#		grep TEST $fname.out >${fname}_test.res
#	done

    echo "---------------------------------------------------------------------"
done