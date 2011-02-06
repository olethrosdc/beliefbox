#data=~/data/regression/maxtemp_one.dat
#b=1;
#time ./bin/kernel_conditional_density_estimation --bandwidth $b --tune_bandwidth --joint --data $data >temp_kde.out;
#grep P_XY temp_kde.out >results/temp_kde.dat;

#time ./bin/kernel_conditional_density_estimation --bandwidth $b --tune_bandwidth --data $data >temp_kcde.out;
#grep P_Y_X temp_kcde.out >results/temp_kcde.dat;





data=~/data/distributions/ring_gauss.dat
b=1
for T in 1000000 #10 100 1000 10000 100000 1000000
do 
    time ./bin/kernel_conditional_density_estimation --bandwidth $b --tune_bandwidth --T $T --joint --data $data >ring_kde_T${T}.out &
    time ./bin/kernel_conditional_density_estimation --bandwidth $b --tune_bandwidth --T $T --data $data >ring_kcde_T${T}.out &
    wait;
    grep P_XY ring_kde_T${T}.out >results/ring_kde_T${T}.dat;
    grep P_Y_X ring_kcde_T${T}.out >results/ring_kcde_T${T}.dat;
done;







