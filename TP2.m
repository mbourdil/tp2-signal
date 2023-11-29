addpath("Data\Data\")
%on teste avec un seul sujet
load('subject_01.mat');
sujet1 = SIGNAL;

%temps de condition 1
debut1 = find(SIGNAL(:,18)==1);
fin1 = debut1+5120;
temps_condition1 = [debut1 fin1];

%temps de condition 2
debut2 = find(SIGNAL(:,19)==1);
fin2 = debut2+5120;
temps_condition2 = [debut2 fin2];

bloc1_cond1 = temps_condition1(1,:);

figure
plot(bloc1_cond1(1):bloc1_cond1(2), sujet1(bloc1_cond1(1):bloc1_cond1(2), 2))
hold on
plot(bloc1_cond1(1):bloc1_cond1(2), sujet1(bloc1_cond1(1):bloc1_cond1(2), 15))
legend('FP1','O1')