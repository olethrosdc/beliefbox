runs=1;
eps=1000;
maze=maze1
steps=100000
for environment in Chain Optimistic ContextBandit;
do
    basedir=$HOME/results/srp
    qsub -v"run=$runs","eps=$eps","environemnt=$environment","maze=$maze","steps=$steps","basedir=$basedir"
done
