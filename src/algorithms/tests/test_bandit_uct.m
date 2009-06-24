## -*- Mode: octave -*-

lw=6;

X=load("results/belief_uct/bandit/uct_is_m.out");
X = sortrows(X,1);
X = sortrows(X,2);
X = sortrows(X,3);
X = sortrows(X,4);
Y=load("results/belief_uct/bandit/baseline_long_double.out");
#Z=load("uct_bnb.out");
gamma = 0.99999;
n_arms=2;
disc = 0;

qa=1;
qb=1;
effort=2;
##effort=9;

hold off;
idx = Y(:,2)==gamma & Y(:,1)==0 & Y(:,3)== n_arms;
r_oracle = Y(idx,6 + disc);

plot([1 2 4 8 16 32], [1 1 1 1 1 1]*(qa*r_oracle - qb*Y(idx,4 + disc)), "0@-;ucb;", "linewidth", lw);

hold on;
idx = Y(:,2)==gamma & Y(:,1)==3 & Y(:,3)== n_arms;
loglog([1 2 4 8 16 32], [1 1 1 1 1 1]*(qa*r_oracle - qb*Y(idx,4 + disc)), "0-;base;", "linewidth", lw);


uct=5 + disc;
oracle=7 + disc;
sel_exp = X(:,3)==gamma & X(:,4)==n_arms;


if (effort==9)
  hold off
endif
r_oracle = mean(X(sel_exp, oracle));
idx = sel_exp & X(:,1)==0; plot(X(idx,effort), qa*r_oracle - qb*X(idx,uct), "1-;serial;", "linewidth", lw);
hold on
idx = sel_exp & X(:,1)==1; plot(X(idx,effort), qa*r_oracle - qb*X(idx,uct), "1@-;random;", "linewidth", lw);
#idx = sel_exp & X(:,1)==2; plot(X(idx,effort), X(idx,oracle)-X(idx,uct), "2-;LB;", "linewidth", lw);
#idx = sel_exp & X(:,1)==; plot(X(idx,effort), qa*r_oracle - qb*X(idx,uct), "2@-;LB {/Symbol g};", "linewidth", lw);
idx = sel_exp & X(:,1)==4; plot(X(idx,effort), qa*r_oracle-qb*X(idx,uct), "2-;thompson;", "linewidth", lw);
#idx = sel_exp & X(:,1)==7; plot(X(idx,effort), qa*r_oracle - qb*X(idx,uct), "3@-;Thompson {/Symbol g};", "linewidth", lw);
#idx = sel_exp & X(:,1)==8; plot(X(idx,effort), X(idx,oracle)-X(idx,uct), "4-;hp UB;", "linewidth", lw);
#idx = sel_exp & X(:,1)==11; plot(X(idx,effort), qa*r_oracle - qb*X(idx,uct), "4@-;x UB {/Symbol g};", "linewidth", lw);
#idx = sel_exp & X(:,1)==13; plot(X(idx,effort), qa*r_oracle - qb*X(idx,uct), "5@-;h.p. UB {/Symbol g};", "linewidth", lw);
idx = sel_exp & X(:,1)==13; plot(X(idx,effort), qa*r_oracle - qb*X(idx,uct), "3-;Alg 2;", "linewidth", lw);
idx = sel_exp & X(:,1)==15; plot(X(idx,effort), qa*r_oracle - qb*X(idx,uct), "3@-;Alg 3;", "linewidth", lw);
idx = sel_exp & X(:,1)==16; plot(X(idx,effort), qa*r_oracle - qb*X(idx,uct), "4-;Alg 3.1;", "linewidth", lw);
idx = sel_exp & X(:,1)==17; plot(X(idx,effort), qa*r_oracle - qb*X(idx,uct), "4@-;Alg 3.2;", "linewidth", lw);


#if (1)
#  sel_exp = Z(:,3)==gamma & Z(:,4)==n_arms;
#  idx = sel_exp & Z(:,1)==14; plot(Z(idx,effort), qa*r_oracle-qb*Z(idx,uct), "3@-;BB;", "linewidth", lw);
#  
#  idx = sel_exp & Z(:,1)==3; plot(Z(idx,effort), qa*r_oracle -qb*Z(idx,uct), "5@-;LB;", "linewidth", lw);
#endif

if (effort==2)
  xlabel("number of expansions");
  axis([1,16]);
else
  xlabel("CPU time");
endif
ylabel("total regret");
##title(strcat("n_{arms}=", n_arms, "{/Symbol g}=", gamma));
legend("location", "southwest");
legend("boxon");
%axis([1,16]);
grid off
##grid("minor")
if (disc) 
  print(strcat("test_bandit_disc_uct_g", num2str(gamma),"_n", num2str(n_arms), ".eps"));
else
  print(strcat("test_bandit_uct_g", num2str(gamma),"_n", num2str(n_arms), ".eps"));
  print(strcat("test_bandit_uct_g", num2str(gamma),"_n", num2str(n_arms), ".fig"));
end
