#! /bin/bash

exe=$HOME/bin/continuous_factored_tree_rl
## this seems to have worked ! ./bin/continuous_factored_tree_rl Pendulum Sarsa 10 10000 16 0.1 0.00 0.5 1 > out
model=BVMM
T=10000
iter=1
erand=0.1
arand=0.0
ulimit -v 2000000
resdir=resdir/half_weights_low_threshold
mkdir -p $resdir

for D in 1 2 4 8 16
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
