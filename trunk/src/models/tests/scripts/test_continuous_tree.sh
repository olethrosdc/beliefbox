#! /bin/bash

exe=$HOME/bin/continuous_factored_tree_rl

model=BVMM
T=1000000
iter=10
erand=0.01
arand=0.01
ulimit -v 2000000
resdir=resdir/half_weights_low_threshold
mkdir -p $resdir

for D in 1 2 3 4 5 6 7 8
do
    for environment in MountainCar MountainCar3D Pendulum LinearBandit
    do
	    echo $T $D $D_c
        cmdline="$exe $environment $model $iter $T $D $erand $arand 0.5 1"
	    echo $cmdline >$resdir/${environment}_d${D}.cmd
	    nice -n 9 $cmdline >$resdir/${environment}_d${D}.out &
    done
    wait;

    for environment in MountainCar MountainCar3D Pendulum LinearBandit
    do
	    grep STATS $resdir/${environment}_d${D}.out >$resdir/${environment}_d${D}_reward.dat
	    grep EPIS $resdir/${environment}_d${D}.out >$resdir/${environment}_d${D}_episode.dat
    done
    wait;
done
