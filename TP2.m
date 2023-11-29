close all
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
bloc1_cond2 = temps_condition2(1,:);

%comparaison signal en occipital et en frontal en condition 1, bloc 1
figure
plot(bloc1_cond1(1):bloc1_cond1(2), sujet1(bloc1_cond1(1):bloc1_cond1(2), 2))
hold on
plot(bloc1_cond1(1):bloc1_cond1(2), sujet1(bloc1_cond1(1):bloc1_cond1(2), 15))
legend('FP1','O1')
title("Comparaison F1 et O1 en condition 1, bloc 1")
hold off

%comparaison signal en occipital et en frontal en condition 2, bloc 1
figure
plot(bloc1_cond2(1):bloc1_cond2(2), sujet1(bloc1_cond2(1):bloc1_cond2(2), 2))
hold on
plot(bloc1_cond2(1):bloc1_cond2(2), sujet1(bloc1_cond2(1):bloc1_cond2(2), 15))
legend('FP1','O1')
title("Comparaison F1 et O1 en condition 2, bloc 1")
hold off


fourier_FP1_cond1 = abs(fft(sujet1(bloc1_cond1(1):bloc1_cond1(2), 2)));
figure
plot(fourier_FP1_cond1(80:120))
res = sum(fourier_FP1_cond1(80:120)) 

% condition 1

% toutes les electrodes O
fourier_O1_bloc1_cond1 = abs(fft(sujet1(bloc1_cond1(1):bloc1_cond1(2), 15)));
fourier_Oz_bloc1_cond1 = abs(fft(sujet1(bloc1_cond1(1):bloc1_cond1(2), 16)));
fourier_O2_bloc1_cond1 = abs(fft(sujet1(bloc1_cond1(1):bloc1_cond1(2), 17)));

% condition 1
for i=1:5
    %identification des blocs
    eval(['bloc',num2str(i),'_cond1 = temps_condition1(',num2str(i),',:);'])
    %calcul des TF
    eval(['fourier_O1_bloc',num2str(i),'_cond1 = abs(fft(sujet1(bloc',num2str(i),'_cond1(1):bloc',num2str(i),'_cond1(2), 15)));'])
    eval(['fourier_Oz_bloc',num2str(i),'_cond1 = abs(fft(sujet1(bloc',num2str(i),'_cond1(1):bloc',num2str(i),'_cond1(2), 16)));'])
    eval(['fourier_O2_bloc',num2str(i),'_cond1 = abs(fft(sujet1(bloc',num2str(i),'_cond1(1):bloc',num2str(i),'_cond1(2), 17)));'])
    %somme pour chaque bloc
    eval(['somme_O1_cond1(',num2str(i),') = sum(fourier_O1_bloc',num2str(i),'_cond1);'])
    eval(['somme_Oz_cond1(',num2str(i),') = sum(fourier_Oz_bloc',num2str(i),'_cond1);'])
    eval(['somme_O2_cond1(',num2str(i),') = sum(fourier_O2_bloc',num2str(i),'_cond1);'])
    %moyenne pour tous les blocs
    mean_O1_cond1 = mean(somme_O1_cond1);
    mean_Oz_cond1 = mean(somme_Oz_cond1);
    mean_O2_cond1 = mean(somme_O2_cond1);
end