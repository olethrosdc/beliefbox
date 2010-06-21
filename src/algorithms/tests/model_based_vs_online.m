figure(1); load model.reward; load sarsa.reward; load qlear.reward;

c=2;
plot(cumsum(sarsa(:,c)), ';sarsa;',cumsum(qlear(:,c)), ';qlear;', cumsum(model(:,c)), ';model;');


figure(2); load model.mse; load sarsa.mse; load qlear.mse;

plot(sarsa(:,1), cumsum(sarsa(:,2)), ';sarsa;',qlear(:,1), cumsum(qlear(:,2)), ';qlear;', model(:,1), cumsum(model(:,2)), ';model;'); 
