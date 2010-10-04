#! /bin/bash

exe=./bin/continuous_factored_tree_rl

environment=Pendulum
model=BVMM
T=1000000
iter=1
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
        cmdline="$exe $environment $model $iter $T $D $erand $arand"
	    echo $cmdline >$resdir/${environment}_d${D}.cmd
	    nice -n 9 $cmdline | grep STATS >$resdir/${environment}_d${D}.out &
    done
    wait;
done
