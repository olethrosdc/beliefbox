ulimit -v 2000000
for i in 1 2 4 8 16 32 64 128 256 512 1024;
do 
    echo $i $j
    time ./bin/conditional_density_estimation --joint --max_depth $i --data ~/data/regression/maxtemp_one.dat >temp_tde_${i}.out; grep P_XY temp_tde_${i}.out >results/temp_tde_${i}.dat;
done;


for i in 1 2 4 8 12
do 
    for j in 1 2 4  8 12
    do 
        echo $i $j
        time ./bin/conditional_density_estimation --n_inputs 1 --max_depth $i --max_depth_cond $j --data ~/data/regression/maxtemp_one.dat >temp_tcde_${i}_${j}.out; 
        grep P_Y_X temp_tcde_${i}_${j}.out >results/temp_tcde_${i}_${j}.dat;
    done;
done




data=~/data/distributions/ring_gauss.dat
b=1
for T in 10 100 1000 10000 100000 1000000
do 
    i=8;
    j=8;
    time ./bin/conditional_density_estimation --n_inputs 1 --max_depth $i --max_depth_cond $j --T $T --data $data --joint >ring_cde_T${T}_${i}_${j}.out; 
    grep P_XY ring_cde_T${T}_${i}_${j}.out >results/ring_cde_T${T}_${i}_${j}.dat;

    time ./bin/conditional_density_estimation --n_inputs 1 --max_depth $i --max_depth_cond $j --T $T --data $data >ring_tcde_T${T}_${i}_${j}.out; 
    grep P_Y_X ring_tcde_T${T}_${i}_${j}.out >results/ring_tcde_T${T}_${i}_${j}.dat;
done;







