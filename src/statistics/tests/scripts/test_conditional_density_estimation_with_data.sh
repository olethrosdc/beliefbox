for i in 1 2 4 8 16 32 64 128 256 512 1024;
do 
    echo $i $j
    time ./bin/conditional_density_estimation --joint --max_depth $i --data ~/data/regression/maxtemp_one.dat >temp_de_${i}.out; grep P_XY temp_de_${i}.out >results/temp_de_${i}.dat;
done;

for i in 1 2 4 8 16 32 64 128 256 512 1024;
do 
    echo $i $j
    time ./bin/conditional_density_estimation --joint --max_depth $i --data ~/data/regression/maxtemp.dat --grid_size 10000 >temp_reg_${i}.out; grep P_XY temp_reg_${i}.out >results/temp_reg_${i}.dat;
done;





for i in 1 2 4 8 16 32;
do 
    for j in 1 2 4 8 16 32;
    do 
        echo $i $j
        time ./bin/conditional_density_estimation --n_inputs 1 --max_depth $i --max_depth_cond $j --data ~/data/regression/maxtemp_one.dat >temp_cde_${i}_${j}.out; grep P_Y_X temp_cde_${i}_${j}.out >results/temp_cde_${i}_${j}.dat;
    done;
done




