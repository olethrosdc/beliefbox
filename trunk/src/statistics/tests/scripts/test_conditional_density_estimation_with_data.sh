
data=~/data/distributions/ring_gauss.dat
b=10
for T in 10 100 1000 10000
do 
    i=8;
    j=8;
    time ./bin/conditional_density_estimation --n_inputs 1 --max_depth $i --max_depth_cond $j --T $T --data $data >ring_tcde_T${T}_${i}_${j}.out; 
    grep P_Y_X ring_tcde_T${T}_${i}_${j}.out >results/ring_tcde_T${T}_${i}_${j}.dat;
    time ./bin/kernel_conditional_density_estimation --bandwidth $b --tune --n_inputs 1 --max_depth $i --max_depth_cond $j --T $T --data $data >ring_kcde_T${T}.out; 
    grep P_Y_X ring_kcde_T${T}.out >results/ring_kcde_T${T}.dat;

    time ./bin/double_kernel_cde --bandwidth $b --tune --n_inputs 1 --max_depth $i --max_depth_cond $j --T $T --data $data >ring_dkcde_T${T}.out; 
    grep P_Y_X ring_dkcde_T${T}.out >results/ring_dkcde_T${T}.dat;
done;



