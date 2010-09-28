datadir=/media/hda2/olethros/data
resdir=~/results/classification
mkdir -p resdir
for dataset in heart ionosphere optdigits pendigits forest spambase usps
do
    train=$datadir/$dataset/valid_train.data
    test=$datadir/$dataset/valid_test.data
    echo ""
    echo "Data: $dataset"
    echo "Training set: $train"
    echo "Test set: $test"
    outdir=$resdir/$dataset/valid/hash_mix
    mkdir -p $outdir
    rm -f $outdir/out
    for (( i=1; i<=32; i++ ))
    do
        ./bin/classify $train $i $test >>$outdir/out
    done
    grep TRAIN $outdir/out >$outdir/train.res
    grep TEST $outdir/out >$outdir/test.res
    echo "---------------------------------------------------------------------"
done