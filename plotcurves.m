clc;
clear;
close all;
%%  Choose the model you want to plot 
% Bound, NDT, Mu, TD
Model = 4;  % 1->TD  2->Bound 3->NDT 4->Mu
ModelNames = {'TD', 'Bound', 'NDT', 'Mu'};
%% Params 
% [-.2, 0] for bound
% [10, 0] for Mu
% [.3, 0] for NDT
% [0.01 , 0] for TD


%%
NumIter = 200;
Cohs = [0:5:25]/100;
AllCoh = repmat(Cohs,NumIter,1);

AllW = [1, 0];  % Choose the Param you want from above section [HCA, LCA]

z = SimModel(Model, AllCoh, AllW, NumIter);
if z
    disp('data was saved successfully')
end
%% Load saved data
load ([ModelNames{Model}, '.mat']);
%% plotting

SmoothBin = 3;
CohBin = Cohs;
LW = 6;
Colors = {[1 0 1], [1 .5 0]};
%% RT
figure, hold on;
for i=1:2
    tmpcolor = Colors{i};
    plot(100*CohBin, mySoomth(AllRT(:, (i)), SmoothBin), 'Color', tmpcolor, 'LineWidth', LW)
    
    
end
set(gcf,'Color','w');
set(gca,'TickDir','out');
set(gca,'FontSize',30);
set(gca,'linewidth',1.5);
set(gca,'Box','off');
title('RT');
ylabel('RT (s)');
xlabel('Stimulus Strength (%coh)')
xlim([0, 26])

%% Conf
figure, hold on;
for i=1:2
    tmpcolor = Colors{i};
    plot(100*CohBin, mySoomth(Allconf(:, (i))/6, SmoothBin), 'Color', tmpcolor, 'LineWidth', LW)
    
    
end
set(gcf,'Color','w');
set(gca,'TickDir','out');
set(gca,'FontSize',30);
set(gca,'linewidth',1.5);
set(gca,'Box','off');
title('Conf');
ylabel('Confidence (a.u.)');
xlabel('Stimulus Strength (%coh)')
xlim([0, 26])
%% ACC
figure, hold on;
for i=1:2
    tmpcolor = Colors{i};
    P{i} = plot(100*CohBin, mySoomth(Decision(:, (i)), SmoothBin), 'Color', tmpcolor, 'LineWidth', LW);
    
    
end
set(gcf,'Color','w');
set(gca,'TickDir','out');
set(gca,'FontSize',30);
set(gca,'linewidth',1.5);
set(gca,'Box','off');
title('Probability Correct');
ylabel('Probability Correct (%)');
xlabel('Stimulus Strength (%coh)')
xlim([0, 26])
% ylim([.35, 1.01])

legend([P{1}, P{2}], {'HCA', 'LCA'})
legend boxoff