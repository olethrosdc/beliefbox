## -*- Mode: octave -*-

lw=2;

X=load("results/belief_uct/bandit/uct_complete.out");
Y=load("results/belief_uct/bandit/baseline_long_double.out");
Z=load("uct_bnb.out");
gamma = 0.9999;
n_arms=2;
disc = 0;

qa=1;
qb=1;

hold off;
idx = Y(:,2)==gamma & Y(:,1)==0 & Y(:,3)== n_arms;
r_oracle = Y(idx,6 + disc);
plot([1 32], [1 1]*(qa*r_oracle - qb*Y(idx,4 + disc)), "0-;ucb;", "linewidth", lw);

hold on;
idx = Y(:,2)==gamma & Y(:,1)==3 & Y(:,3)== n_arms;
plot([1 32], [1 1]*(qa*r_oracle - qb*Y(idx,4 + disc)), "0@-;base;", "linewidth", lw);


uct=5 + disc;
oracle=7 + disc;
sel_exp = X(:,3)==gamma & X(:,4)==n_arms;


idx = sel_exp & X(:,1)==0; plot(X(idx,2), qa*r_oracle - qb*X(idx,uct), "1-;serial;", "linewidth", lw);
idx = sel_exp & X(:,1)==1; plot(X(idx,2), qa*r_oracle - qb*X(idx,uct), "1@-;random;", "linewidth", lw);
#idx = sel_exp & X(:,1)==2; plot(X(idx,2), X(idx,oracle)-X(idx,uct), "2-;LB;", "linewidth", lw);
idx = sel_exp & X(:,1)==3; plot(X(idx,2), qa*r_oracle - qb*X(idx,uct), "2@-;LB {/Symbol g};", "linewidth", lw);
%idx = sel_exp & X(:,1)==6; plot(X(idx,2), X(idx,oracle)-X(idx,uct), "3-;thompson;", "linewidth", lw);
idx = sel_exp & X(:,1)==7; plot(X(idx,2), qa*r_oracle - qb*X(idx,uct), "3@-;Thompson {/Symbol g};", "linewidth", lw);
#idx = sel_exp & X(:,1)==8; plot(X(idx,2), X(idx,oracle)-X(idx,uct), "4-;hp UB;", "linewidth", lw);
idx = sel_exp & X(:,1)==11; plot(X(idx,2), qa*r_oracle - qb*X(idx,uct), "4@-;x UB {/Symbol g};", "linewidth", lw);
idx = sel_exp & X(:,1)==13; plot(X(idx,2), qa*r_oracle - qb*X(idx,uct), "5@-;h.p. UB {/Symbol g};", "linewidth", lw);


if (1)
  sel_exp = Z(:,3)==gamma & Z(:,4)==n_arms;
  idx = sel_exp & Z(:,1)==14; plot(Z(idx,2), qa*r_oracle-qb*Z(idx,uct), "3@-;BB;", "linewidth", lw);
  
  idx = sel_exp & Z(:,1)==3; plot(Z(idx,2), qa*r_oracle -qb*Z(idx,uct), "5@-;LB;", "linewidth", lw);
endif

xlabel("number of expansions");
ylabel("total regret");
##title(strcat("n_{arms}=", n_arms, "{/Symbol g}=", gamma));
legend("location", "southwest");
legend("boxon");
axis([1,16]);
grid on
#grid("minor")
#if (disc) 
#  print(strcat("test_bandit_disc_uct_g", num2str(gamma),"_n", num2str(n_arms), ".eps"));
#else
#  print(strcat("test_bandit_uct_g", num2str(gamma),"_n", num2str(n_arms), ".eps"));
#end
