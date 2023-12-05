%% Pré traitement des données 
close all
addpath("Data\Data\")

%on s'interesse au sujet 1
load('subject_01.mat');

%on détermine les temps où il est en condition 1 ou 2

%temps de condition 1
debut1 = find(SIGNAL(:,18)==1);
fin1 = debut1+5120;

%temps de condition 2
debut2 = find(SIGNAL(:,19)==1);
fin2 = debut2+5120;

%% Comparaison du signal entre occipital et frontal 


bloc1_cond1 = [debut1, fin1];
bloc1_cond2 = [debut2, fin2];

% sujet 1, condition 1, bloc 1
figure
plot(bloc1_cond1(1):bloc1_cond1(2), sujet1(bloc1_cond1(1):bloc1_cond1(2), 2))
hold on
plot(bloc1_cond1(1):bloc1_cond1(2), sujet1(bloc1_cond1(1):bloc1_cond1(2), 15))
legend('FP1','O1')
title("Comparaison F1 et O1 pour le sujet 1 en condition 1, bloc 1")
hold off

%sujet 1, condition 2, bloc 1
figure
plot(bloc1_cond2(1):bloc1_cond2(2), sujet1(bloc1_cond2(1):bloc1_cond2(2), 2))
hold on
plot(bloc1_cond2(1):bloc1_cond2(2), sujet1(bloc1_cond2(1):bloc1_cond2(2), 15))
legend('FP1','O1')
title("Comparaison F1 et O1 pour le sujet 1 en condition 2, bloc 1")
hold off

%le signal prend bien des valeurs plus importantes en région occipitale

%% TF d'un bloc

fourier_O1_cond1 = abs(fft(sujet1(bloc1_cond1(1):bloc1_cond1(2), 15)));
figure
plot(80:120, fourier_O1_cond1(80:120))
title("Transformée de Fourier de l'électrode O1, sujet 1, condition 1, bloc 1")

%% Comparaison des électrodes O en condition 1 et condition 2 pour le sujet 1

% condition 1
for i=1:5
    %identification des blocs
    eval(['bloc',num2str(i),'_cond1 = temps_condition1(',num2str(i),',:);'])

    %normalisation du signal
    eval(['O1_bloc',num2str(i),'_cond1 = sujet1(bloc',num2str(i),'_cond1(1):bloc',num2str(i),'_cond1(2), 15);'])
    eval(['norm_O1_bloc',num2str(i),'__cond1 = ((O1_bloc',num2str(i),'_cond1) - min(O1_bloc',num2str(i),'_cond1)) / (max(O1_bloc',num2str(i),'_cond1)-min(O1_bloc',num2str(i),'_cond1));'])
    eval(['Oz_bloc',num2str(i),'_cond1 = sujet1(bloc',num2str(i),'_cond1(1):bloc',num2str(i),'_cond1(2), 16);'])
    eval(['norm_Oz_bloc',num2str(i),'__cond1 = ((Oz_bloc',num2str(i),'_cond1) - min(Oz_bloc',num2str(i),'_cond1)) / (max(Oz_bloc',num2str(i),'_cond1)-min(Oz_bloc',num2str(i),'_cond1));'])
    eval(['O2_bloc',num2str(i),'_cond1 = sujet1(bloc',num2str(i),'_cond1(1):bloc',num2str(i),'_cond1(2), 17);'])
    eval(['norm_O2_bloc',num2str(i),'__cond1 = ((O2_bloc',num2str(i),'_cond1) - min(O2_bloc',num2str(i),'_cond1)) / (max(O2_bloc',num2str(i),'_cond1)-min(O2_bloc',num2str(i),'_cond1));'])
    
    %calcul des TF
    eval(['fourier_O1_bloc',num2str(i),'_cond1 = abs(fft(norm_O1_bloc',num2str(i),'__cond1));'])
    eval(['fourier_Oz_bloc',num2str(i),'_cond1 = abs(fft(norm_Oz_bloc',num2str(i),'__cond1));'])
    eval(['fourier_O2_bloc',num2str(i),'_cond1 = abs(fft(norm_O2_bloc',num2str(i),'__cond1));'])
    
    %somme pour chaque bloc
    eval(['somme_O1_cond1(',num2str(i),') = sum(fourier_O1_bloc',num2str(i),'_cond1);'])
    eval(['somme_Oz_cond1(',num2str(i),') = sum(fourier_Oz_bloc',num2str(i),'_cond1);'])
    eval(['somme_O2_cond1(',num2str(i),') = sum(fourier_O2_bloc',num2str(i),'_cond1);'])
    
    %moyenne pour tous les blocs
    mean_O1_cond1 = mean(somme_O1_cond1);
    mean_Oz_cond1 = mean(somme_Oz_cond1);
    mean_O2_cond1 = mean(somme_O2_cond1);
    mean_O_cond1 = [mean_O1_cond1 mean_Oz_cond1 mean_O2_cond1];
end

disp('Moyenne pour les electrodes O sur le sujet 1 en condition 1 = ')
disp(mean(mean_O_cond1))

%condition 2

for i=1:5
    %identification des blocs
    eval(['bloc',num2str(i),'_cond2 = temps_condition2(',num2str(i),',:);'])

    %normalisation du signal
    eval(['O1_bloc',num2str(i),'_cond2 = sujet1(bloc',num2str(i),'_cond2(1):bloc',num2str(i),'_cond2(2), 15);'])
    eval(['norm_O1_bloc',num2str(i),'__cond2 = ((O1_bloc',num2str(i),'_cond2) - min(O1_bloc',num2str(i),'_cond2)) / (max(O1_bloc',num2str(i),'_cond2)-min(O1_bloc',num2str(i),'_cond2));'])
    eval(['Oz_bloc',num2str(i),'_cond2 = sujet1(bloc',num2str(i),'_cond2(1):bloc',num2str(i),'_cond2(2), 16);'])
    eval(['norm_Oz_bloc',num2str(i),'__cond2 = ((Oz_bloc',num2str(i),'_cond2) - min(Oz_bloc',num2str(i),'_cond2)) / (max(Oz_bloc',num2str(i),'_cond2)-min(Oz_bloc',num2str(i),'_cond2));'])
    eval(['O2_bloc',num2str(i),'_cond2 = sujet1(bloc',num2str(i),'_cond2(1):bloc',num2str(i),'_cond2(2), 17);'])
    eval(['norm_O2_bloc',num2str(i),'__cond2 = ((O2_bloc',num2str(i),'_cond2) - min(O2_bloc',num2str(i),'_cond2)) / (max(O2_bloc',num2str(i),'_cond2)-min(O2_bloc',num2str(i),'_cond2));'])
    
    %calcul des TF
    eval(['fourier_O1_bloc',num2str(i),'_cond2 = abs(fft(norm_O1_bloc',num2str(i),'__cond2));'])
    eval(['fourier_Oz_bloc',num2str(i),'_cond2 = abs(fft(norm_Oz_bloc',num2str(i),'__cond2));'])
    eval(['fourier_O2_bloc',num2str(i),'_cond2 = abs(fft(norm_O2_bloc',num2str(i),'__cond2));'])
    
    %somme pour chaque bloc
    eval(['somme_O1_cond2(',num2str(i),') = sum(fourier_O1_bloc',num2str(i),'_cond2);'])
    eval(['somme_Oz_cond2(',num2str(i),') = sum(fourier_Oz_bloc',num2str(i),'_cond2);'])
    eval(['somme_O2_cond2(',num2str(i),') = sum(fourier_O2_bloc',num2str(i),'_cond2);'])
    
    %moyenne pour tous les blocs
    mean_O1_cond2 = mean(somme_O1_cond2);
    mean_Oz_cond2 = mean(somme_Oz_cond2);
    mean_O2_cond2 = mean(somme_O2_cond2);
    mean_O_cond2 = [mean_O1_cond2 mean_Oz_cond2 mean_O2_cond2];
end

disp('Moyenne pour les electrodes O sur le sujet 1 en condition 2 = ')
disp(mean(mean_O_cond2))

%en moyenne plus fort en condition 1 que 2

%% Maintenant avec tous les sujets et toutes les electrodes

% initialisation des matrices de stockage des moyennes des électrodes O 
% selon les sujets en ligne et les conditions en colonne.

moyennes_O = zeros(20,2);
moyennes_F = zeros(20,2);
moyennes_T = zeros(20,2);
moyennes_P = zeros(20,2);

for k=1:20
    if k<10
        eval(['load(''subject_0',num2str(k),'.mat'');'])
    end

    if k>=10
        eval(['load(''subject_',num2str(k),'.mat'');'])
    end

    %temps de condition 1
    debut1 = find(SIGNAL(:,18)==1);
    fin1 = debut1+5120;
    temps_condition1 = [debut1 fin1];
    
    %temps de condition 2
    debut2 = find(SIGNAL(:,19)==1);
    fin2 = debut2+5120;
    temps_condition2 = [debut2 fin2];


% condition 1
% pour tous les blocs
for i=1:5
    %identification des blocs
    eval(['bloc',num2str(i),'_cond1 = temps_condition1(',num2str(i),',:);'])

    %normalisation du signal

    eval(['O1_bloc',num2str(i),'_cond1 = SIGNAL(bloc',num2str(i),'_cond1(1):bloc',num2str(i),'_cond1(2), 15);'])
    eval(['norm_O1_bloc',num2str(i),'__cond1 = ((O1_bloc',num2str(i),'_cond1) - min(O1_bloc',num2str(i),'_cond1)) / (max(O1_bloc',num2str(i),'_cond1)-min(O1_bloc',num2str(i),'_cond1));'])
    eval(['O1_bloc',num2str(i),'_cond1 = SIGNAL(bloc',num2str(i),'_cond1(1):bloc',num2str(i),'_cond1(2), 15);'])
    eval(['norm_O1_bloc',num2str(i),'__cond1 = ((O1_bloc',num2str(i),'_cond1) - min(O1_bloc',num2str(i),'_cond1)) / (max(O1_bloc',num2str(i),'_cond1)-min(O1_bloc',num2str(i),'_cond1));'])
    eval(['O1_bloc',num2str(i),'_cond1 = SIGNAL(bloc',num2str(i),'_cond1(1):bloc',num2str(i),'_cond1(2), 15);'])
    eval(['norm_O1_bloc',num2str(i),'__cond1 = ((O1_bloc',num2str(i),'_cond1) - min(O1_bloc',num2str(i),'_cond1)) / (max(O1_bloc',num2str(i),'_cond1)-min(O1_bloc',num2str(i),'_cond1));'])
    eval(['O1_bloc',num2str(i),'_cond1 = SIGNAL(bloc',num2str(i),'_cond1(1):bloc',num2str(i),'_cond1(2), 15);'])
    eval(['norm_O1_bloc',num2str(i),'__cond1 = ((O1_bloc',num2str(i),'_cond1) - min(O1_bloc',num2str(i),'_cond1)) / (max(O1_bloc',num2str(i),'_cond1)-min(O1_bloc',num2str(i),'_cond1));'])
    eval(['O1_bloc',num2str(i),'_cond1 = SIGNAL(bloc',num2str(i),'_cond1(1):bloc',num2str(i),'_cond1(2), 15);'])
    eval(['norm_O1_bloc',num2str(i),'__cond1 = ((O1_bloc',num2str(i),'_cond1) - min(O1_bloc',num2str(i),'_cond1)) / (max(O1_bloc',num2str(i),'_cond1)-min(O1_bloc',num2str(i),'_cond1));'])
    eval(['O1_bloc',num2str(i),'_cond1 = SIGNAL(bloc',num2str(i),'_cond1(1):bloc',num2str(i),'_cond1(2), 15);'])
    eval(['norm_O1_bloc',num2str(i),'__cond1 = ((O1_bloc',num2str(i),'_cond1) - min(O1_bloc',num2str(i),'_cond1)) / (max(O1_bloc',num2str(i),'_cond1)-min(O1_bloc',num2str(i),'_cond1));'])
    eval(['O1_bloc',num2str(i),'_cond1 = SIGNAL(bloc',num2str(i),'_cond1(1):bloc',num2str(i),'_cond1(2), 15);'])
    eval(['norm_O1_bloc',num2str(i),'__cond1 = ((O1_bloc',num2str(i),'_cond1) - min(O1_bloc',num2str(i),'_cond1)) / (max(O1_bloc',num2str(i),'_cond1)-min(O1_bloc',num2str(i),'_cond1));'])
    
    eval(['O1_bloc',num2str(i),'_cond1 = SIGNAL(bloc',num2str(i),'_cond1(1):bloc',num2str(i),'_cond1(2), 15);'])
    eval(['norm_O1_bloc',num2str(i),'__cond1 = ((O1_bloc',num2str(i),'_cond1) - min(O1_bloc',num2str(i),'_cond1)) / (max(O1_bloc',num2str(i),'_cond1)-min(O1_bloc',num2str(i),'_cond1));'])
    eval(['Oz_bloc',num2str(i),'_cond1 = SIGNAL(bloc',num2str(i),'_cond1(1):bloc',num2str(i),'_cond1(2), 16);'])
    eval(['norm_Oz_bloc',num2str(i),'__cond1 = ((Oz_bloc',num2str(i),'_cond1) - min(Oz_bloc',num2str(i),'_cond1)) / (max(Oz_bloc',num2str(i),'_cond1)-min(Oz_bloc',num2str(i),'_cond1));'])
    eval(['O2_bloc',num2str(i),'_cond1 = SIGNAL(bloc',num2str(i),'_cond1(1):bloc',num2str(i),'_cond1(2), 17);'])
    eval(['norm_O2_bloc',num2str(i),'__cond1 = ((O2_bloc',num2str(i),'_cond1) - min(O2_bloc',num2str(i),'_cond1)) / (max(O2_bloc',num2str(i),'_cond1)-min(O2_bloc',num2str(i),'_cond1));'])
    
    %calcul des TF
    eval(['fourier_O1_bloc',num2str(i),'_cond1 = abs(fft(norm_O1_bloc',num2str(i),'__cond1));'])
    eval(['fourier_Oz_bloc',num2str(i),'_cond1 = abs(fft(norm_Oz_bloc',num2str(i),'__cond1));'])
    eval(['fourier_O2_bloc',num2str(i),'_cond1 = abs(fft(norm_O2_bloc',num2str(i),'__cond1));'])

    eval(['fourier_O1_bloc',num2str(i),'_cond1 = fourier_O1_bloc',num2str(i),'_cond1(80:120);'])
    eval(['fourier_Oz_bloc',num2str(i),'_cond1 = fourier_Oz_bloc',num2str(i),'_cond1(80:120);'])
    eval(['fourier_O2_bloc',num2str(i),'_cond1 = fourier_O2_bloc',num2str(i),'_cond1(80:120);'])
    
    %somme pour chaque bloc
    eval(['somme_O1_cond1(',num2str(i),') = sum(fourier_O1_bloc',num2str(i),'_cond1);'])
    eval(['somme_Oz_cond1(',num2str(i),') = sum(fourier_Oz_bloc',num2str(i),'_cond1);'])
    eval(['somme_O2_cond1(',num2str(i),') = sum(fourier_O2_bloc',num2str(i),'_cond1);'])
    
    %moyenne pour tous les blocs
    mean_O1_cond1 = mean(somme_O1_cond1);
    mean_Oz_cond1 = mean(somme_Oz_cond1);
    mean_O2_cond1 = mean(somme_O2_cond1);
    mean_O_cond1 = [mean_O1_cond1 mean_Oz_cond1 mean_O2_cond1];
end

moyennes_O(k,1) = mean(mean_O_cond1);

%condition 2

for i=1:5
    %identification des blocs
    eval(['bloc',num2str(i),'_cond2 = temps_condition2(',num2str(i),',:);'])

    %normalisation du signal
    eval(['O1_bloc',num2str(i),'_cond2 = SIGNAL(bloc',num2str(i),'_cond2(1):bloc',num2str(i),'_cond2(2), 15);'])
    eval(['norm_O1_bloc',num2str(i),'__cond2 = ((O1_bloc',num2str(i),'_cond2) - min(O1_bloc',num2str(i),'_cond2)) / (max(O1_bloc',num2str(i),'_cond2)-min(O1_bloc',num2str(i),'_cond2));'])
    eval(['Oz_bloc',num2str(i),'_cond2 = SIGNAL(bloc',num2str(i),'_cond2(1):bloc',num2str(i),'_cond2(2), 16);'])
    eval(['norm_Oz_bloc',num2str(i),'__cond2 = ((Oz_bloc',num2str(i),'_cond2) - min(Oz_bloc',num2str(i),'_cond2)) / (max(Oz_bloc',num2str(i),'_cond2)-min(Oz_bloc',num2str(i),'_cond2));'])
    eval(['O2_bloc',num2str(i),'_cond2 = SIGNAL(bloc',num2str(i),'_cond2(1):bloc',num2str(i),'_cond2(2), 17);'])
    eval(['norm_O2_bloc',num2str(i),'__cond2 = ((O2_bloc',num2str(i),'_cond2) - min(O2_bloc',num2str(i),'_cond2)) / (max(O2_bloc',num2str(i),'_cond2)-min(O2_bloc',num2str(i),'_cond2));'])
    
    %calcul des TF
    eval(['fourier_O1_bloc',num2str(i),'_cond2 = abs(fft(norm_O1_bloc',num2str(i),'__cond2));'])
    eval(['fourier_Oz_bloc',num2str(i),'_cond2 = abs(fft(norm_Oz_bloc',num2str(i),'__cond2));'])
    eval(['fourier_O2_bloc',num2str(i),'_cond2 = abs(fft(norm_O2_bloc',num2str(i),'__cond2));'])
    
    %on ne prend que les fréquences entre 80 et 120 
    eval(['fourier_O1_bloc',num2str(i),'_cond2 = fourier_O1_bloc',num2str(i),'_cond2(80:120);'])
    eval(['fourier_Oz_bloc',num2str(i),'_cond2 = fourier_Oz_bloc',num2str(i),'_cond2(80:120);'])
    eval(['fourier_O2_bloc',num2str(i),'_cond2 = fourier_O2_bloc',num2str(i),'_cond2(80:120);'])



    %somme pour chaque bloc
    eval(['somme_O1_cond2(',num2str(i),') = sum(fourier_O1_bloc',num2str(i),'_cond2);'])
    eval(['somme_Oz_cond2(',num2str(i),') = sum(fourier_Oz_bloc',num2str(i),'_cond2);'])
    eval(['somme_O2_cond2(',num2str(i),') = sum(fourier_O2_bloc',num2str(i),'_cond2);'])
    
    %moyenne pour tous les blocs
    mean_O1_cond2 = mean(somme_O1_cond2);
    mean_Oz_cond2 = mean(somme_Oz_cond2);
    mean_O2_cond2 = mean(somme_O2_cond2);
    mean_O_cond2 = [mean_O1_cond2 mean_Oz_cond2 mean_O2_cond2];
end

moyennes_O(k,2) = mean(mean_O_cond2);

end

%% Représentation


%Toutes les O
figure
plot(1:1:20,moyennes_O(:,1),'*')
hold on
plot(1:1:20,moyennes_O(:,2),'o')
legend('condition 1 = fermés','condition 2 = ouverts')
title('moyennes des électrodes O selon les sujets et les conditions')
xlabel('sujets')
ylabel('Hz')


