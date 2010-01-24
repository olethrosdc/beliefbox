T=10000;

fname=out1
rm -f $fname
./bin/bayesian_markov_chain_text data/alice.txt 0 $T 1 >>$fname
./bin/bayesian_markov_chain_text data/alice.txt 1 $T 0.5 >>$fname
./bin/bayesian_markov_chain_text data/alice.txt 2 $T 0.292893 >>$fname
./bin/bayesian_markov_chain_text data/alice.txt 3 $T 0.206299 >>$fname
./bin/bayesian_markov_chain_text data/alice.txt 4 $T 0.159104 >>$fname
./bin/bayesian_markov_chain_text data/alice.txt 5 $T 0.129449 >>$fname

fname=out2
rm -f $fname
./bin/bayesian_markov_chain_text data/alice.txt 0 $T 0.5 >>$fname
./bin/bayesian_markov_chain_text data/alice.txt 1 $T 0.292893 >>$fname
./bin/bayesian_markov_chain_text data/alice.txt 2 $T 0.206299 >>$fname
./bin/bayesian_markov_chain_text data/alice.txt 3 $T 0.159104 >>$fname
./bin/bayesian_markov_chain_text data/alice.txt 4 $T 0.129449 >>$fname
./bin/bayesian_markov_chain_text data/alice.txt 5 $T 0.109101 >>out

fname=out3
rm -f $fname
./bin/bayesian_markov_chain_text data/alice.txt 0 $T 0.5 >>$fname
./bin/bayesian_markov_chain_text data/alice.txt 1 $T 0.5 >>$fname
./bin/bayesian_markov_chain_text data/alice.txt 2 $T 0.5 >>$fname
./bin/bayesian_markov_chain_text data/alice.txt 3 $T 0.5 >>$fname
./bin/bayesian_markov_chain_text data/alice.txt 4 $T 0.5 >>$fname
./bin/bayesian_markov_chain_text data/alice.txt 5 $T 0.5 >>out

fname=out4
rm -f $fname
./bin/bayesian_markov_chain_text data/alice.txt 0 $T 0.027 >>$fname
./bin/bayesian_markov_chain_text data/alice.txt 1 $T 0.027 >>$fname
./bin/bayesian_markov_chain_text data/alice.txt 2 $T 0.027 >>$fname
./bin/bayesian_markov_chain_text data/alice.txt 3 $T 0.027 >>$fname
./bin/bayesian_markov_chain_text data/alice.txt 4 $T 0.027 >>$fname
./bin/bayesian_markov_chain_text data/alice.txt 5 $T 0.027 >>out

fname=out5
rm -f $fname
./bin/bayesian_markov_chain_text data/alice.txt 0 $T 1 >>$fname
./bin/bayesian_markov_chain_text data/alice.txt 1 $T 0.027 >>$fname
./bin/bayesian_markov_chain_text data/alice.txt 2 $T 0.0013594 >>$fname
./bin/bayesian_markov_chain_text data/alice.txt 3 $T 0.0090822 >>$fname
./bin/bayesian_markov_chain_text data/alice.txt 4 $T 0.0068194 >>$fname
./bin/bayesian_markov_chain_text data/alice.txt 5 $T 0.0054593 >>out

## 0.500000   0.292893   0.206299   0.159104   0.129449   0.109101   0.094276
