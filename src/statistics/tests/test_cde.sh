for d in 1 2 4 8 16 32
do
    for dc in 1 2 4 8 16 32
    do
        ./bin/conditional_density_estimation  --data ~/data/regression/maxtemp_one.dat --max_depth $d --max_depth_cond $dc | grep P_Y_X > P_Y_X_${d}_${dc}
    done
done