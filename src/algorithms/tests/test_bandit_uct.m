## -*- Mode: octave -*-

lw=2;

X=load("results/belief_uct/bandit/uct_complete.out");
Y=load("results/belief_uct/bandit/baseline_long_double.out");
Z=load("uct_bnb.out");
gamma = 0.999;
n_arms=2;
disc = 0;

hold off;
idx = Y(:,2)==gamma & Y(:,1)==0 & Y(:,3)== n_arms;
plot([1 32], [1 1]*(Y(idx, 6 + disc) - Y(idx,4 + disc)), "0-;ucb;", "linewidth", lw);

hold on;
idx = Y(:,2)==gamma & Y(:,1)==3 & Y(:,3)== n_arms;
plot([1 32], [1 1]*(Y(idx, 6 + disc) - Y(idx,4 + disc)), "0@-;base;", "linewidth", lw);


uct=5 + disc;
oracle=7 + disc;
sel_exp = X(:,3)==gamma & X(:,4)==n_arms;

idx = sel_exp & X(:,1)==0; plot(X(idx,2), X(idx,oracle)-X(idx,uct), "1-;serial;", "linewidth", lw);
idx = sel_exp & X(:,1)==1; plot(X(idx,2), X(idx,oracle)-X(idx,uct), "1@-;random;", "linewidth", lw);
idx = sel_exp & X(:,1)==2; plot(X(idx,2), X(idx,oracle)-X(idx,uct), "2-;LB;", "linewidth", lw);
idx = sel_exp & X(:,1)==3; plot(X(idx,2), X(idx,oracle)-X(idx,uct), "2@-;LB gamma;", "linewidth", lw);
#idx = sel_exp & X(:,1)==6; plot(X(idx,2), X(idx,oracle)-X(idx,uct), "3-;thompson;", "linewidth", lw);
#idx = sel_exp & X(:,1)==7; plot(X(idx,2), X(idx,oracle)-X(idx,uct), "3@-;thompson gamma;", "linewidth", lw);
idx = sel_exp & X(:,1)==8; plot(X(idx,2), X(idx,oracle)-X(idx,uct), "4-;hp UB;", "linewidth", lw);
idx = sel_exp & X(:,1)==9; plot(X(idx,2), X(idx,oracle)-X(idx,uct), "4@-;hp UB gamma;", "linewidth", lw);

sel_exp = Z(:,3)==gamma & Z(:,4)==n_arms;


idx = sel_exp & Z(:,1)==14; semilogy(Z(idx,2), Z(idx,oracle)-Z(idx,uct), "3@-;BB;", "linewidth", lw);

idx = sel_exp & Z(:,1)==3; semilogy(Z(idx,2), Z(idx,oracle)-Z(idx,uct), "5@-;LB;", "linewidth", lw);

xlabel("number of expansions");
ylabel("total regret");
title(strcat("n_{arms}=", n_arms, "{\Sym gamma}=", gamma));
legend("location", "northeast");
legend("boxon");
grid on
print("test_bandit_uct.eps");
