datadir=~/data/classification
resdir=~/results/classification
mkdir -p resdir
for dataset in optdigits pendigits spect
do
    train=$datadir/$dataset/train.dat
    test=$datadir/$dataset/test.dat
    echo ""
    echo "Data: $dataset"
    echo "Training set: $train"
    echo "Test set: $test"
    outdir=$resdir/$dataset/hash_mix
    mkdir -p $outdir
	for knn in 1 2 3 4
	do
		fname=$outdir/knn${knn}
		/usr/bin/time --verbose --output $fname.time ./bin/classify --train $train --test $test --randomise --normalise --knn $knn  >$fname.out
		grep TRAIN $fname.out >${fname}_train.res
		grep TEST $fname.out >${fname}_test.res
		for mix in 1 2 4 8 16
		do
			fname=$outdir/knn${knn}_mix${mix}
			/usr/bin/time --verbose --output $fname.time ./bin/classify --train $train --test $test --randomise --normalise --knn $knn --mixture $mix >$fname.out
			grep TRAIN $fname.out >${fname}_train.res
			grep TEST $fname.out >${fname}_test.res
		done
	done

	fname=$outdir/tree${knn}
	for tree in 1 2 4 8 16 32 64 128
	do
		fname=$outdir/tree${tree}
		/usr/bin/time --verbose --output $fname.time ./bin/classify --train $train --test $test --randomise --normalise --tree $tree >$fname.out
		grep TRAIN $fname.out >${fname}_train.res
		grep TEST $fname.out >${fname}_test.res
	done
    echo "---------------------------------------------------------------------"
done