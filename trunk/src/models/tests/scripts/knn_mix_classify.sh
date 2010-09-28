datadir=/media/hda2/olethros/data
resdir=~/results/classification
mkdir -p resdir
for dataset in usps #heart ionosphere optdigits pendigits forest spambase usps
do
    train=$datadir/$dataset/valid_train.data
    test=$datadir/$dataset/valid_test.data
    echo ""
    echo "Data: $dataset"
    echo "Training set: $train"
    echo "Test set: $test"
    for (( k=1; k<=8; k++ ))
    do
        outdir=$resdir/$dataset/valid/knn${k}
        mkdir -p $outdir
        rm -f $outdir/*
        for (( i=1; i<=32; i++ ))
        do
            ./bin/knn_mix_classify $train $k $i $test >>$outdir/out
        done
        grep TRAIN $outdir/out >$outdir/train.res
        grep TEST $outdir/out >$outdir/test.res
    done
    echo "---------------------------------------------------------------------"
done
