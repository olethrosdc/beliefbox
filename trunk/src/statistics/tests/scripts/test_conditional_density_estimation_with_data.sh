name=gaussian_mix
data=~/data/distributions/${name}_train.dat
test_data=~/data/distributions/${name}_test.dat
n_inputs=1
b=10
for T in 10 20 50 100 200 500 1000 2000 5000 10000 20000 50000 100000
do 
    i=8;
    j=8;
    time ./bin/conditional_density_estimation --n_inputs $n_inputs --max_depth $i --max_depth_cond $j --T $T --data $data --test $test_data >${name}_tcde_T${T}_${i}_${j}.out; 
    grep P_Y_X ${name}_tcde_T${T}_${i}_${j}.out >results/${name}_tcde_T${T}_${i}_${j}.dat;
    time ./bin/kernel_conditional_density_estimation --bandwidth $b --tune --n_inputs $n_inputs --max_depth $i --max_depth_cond $j --T $T --data $data --test $test_data >${name}_kcde_T${T}.out; 
    grep P_Y_X ${name}_kcde_T${T}.out >results/${name}_kcde_T${T}.dat;

    time ./bin/double_kernel_cde --bandwidth $b --tune --n_inputs $n_inputs --max_depth $i --max_depth_cond $j --T $T --data $data --test $test_data >${name}_dkcde_T${T}.out; 
    grep P_Y_X ${name}_dkcde_T${T}.out >results/${name}_dkcde_T${T}.dat;
done;



