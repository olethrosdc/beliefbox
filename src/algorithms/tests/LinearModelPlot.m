clc
clear all
close all

actions = 3;

input_target_file = 'Input_samples.txt';
X = textread(input_target_file);
X = X(:,1:end-1);
hold on;

output_target_file = ['Output_samples_action_' num2str(0)];
Y    = textread(output_target_file);
Y    = Y(:,1:end-1);
plot(X(:,1), Y(:,1),'*','MarkerSize',7);
hold on;

marginal_target_file = 'Predicted_Output_Marginal_Sampling';
Ymarginal = textread(marginal_target_file);
Ymarginal   = Ymarginal(:,1:end-1);
plot(X(:,1), Ymarginal(:,1),'*','Color','r','MarkerSize',7);
hold on;

tompson_target_file_1 = 'Predicted_Output_Thompson_Sampling_1';
Ytompson1 = textread(tompson_target_file_1);
Ytompson1   = Ytompson1(:,1:end-1);
plot(X(:,1), Ytompson1(:,1),'*','Color','y','MarkerSize',7);
hold on
tompson_target_file_2 = 'Predicted_Output_Thompson_Sampling_2';
Ytompson2 = textread(tompson_target_file_2);
Ytompson2   = Ytompson2(:,1:end-1);
plot(X(:,1), Ytompson2(:,1),'*','Color','g','MarkerSize',7);

legend('Expected','Marginal','Sampling 1','Sampling 2');