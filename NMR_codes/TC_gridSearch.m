% Run Timur-Coates Model Estimates


% Range over pairs of m and n values
close all
clear

load enso 

%siteList = [{'Site1-WellG6above'} {'Site1-WellG6below'} {'Site1-WellG5above'}...
%    {'Site1-WellG5below'} {'Site2-WellPN1'} {'Site2-WellPN2'} {'Site1-WellG6'} {'Site1-WellG5'}];

siteList = [{'dpnmr_larned_east'} {'dpnmr_larned_lwph'}]% {'dpnmr_larned_west'}...
    %{'dpnmr_leque_east'},{'dpnmr_leque_west'}];

%siteList = [{'dpnmrA11'} {'dpnmrA12'} {'dpnmrC1S'} {'dpnmrC1SE'} {'dpnmrC1SW'}];

% Matrix Structure... Repeat for each site
%       Bootstrap           
% c
% m
% n

%cutoff = [20 33 40 50 60 80 100 120 140 160 200 240 280 320 380 440 500 560 620]*10^-3;
cutoff = (20:2:800)*10^-3;
%cutoff = 510*10^-3;
m = [];
n = [];

figureson = 0;
wDirect = 1;
currentRow = -3;
%%
tic
totalErrorMatrix = zeros(length(siteList),length(cutoff));
totalmMatrix = zeros(length(siteList),length(cutoff));
totalnMatrix = zeros(length(siteList),length(cutoff));
totalcMatrix = zeros(length(siteList),length(cutoff));

for j = 1:length(siteList)
    currentFitMatrix = [];
    currentErrorMatrix = [];
    baseName = siteList{j};
 
    mTemp = zeros(1,length(cutoff));
    nTemp = zeros(1,length(cutoff));
    cTemp = zeros(1,length(cutoff));
    errorTemp = zeros(1,length(cutoff));
    
    
    for i = 1:length(cutoff)
 
        [K,z,T2dist,T2logbins,kTC_best,bestFitMatrix,totalErrorEstimate] = computeTCperm_2(baseName,n,m,cutoff(i),figureson);
 
        mTemp(i) = bestFitMatrix(3);
        nTemp(i) = bestFitMatrix(2);
        cTemp(i) = bestFitMatrix(1);
        errorTemp(i) = totalErrorEstimate;
    end
    
    totalmMatrix(j,:) = mTemp;
    totalnMatrix(j,:) = nTemp;
    totalcMatrix(j,:) = cTemp;
    totalErrorMatrix(j,:) = errorTemp;
        
end
toc
 
save('optimalCutoffTableKansas_var_n_m.mat','totalmMatrix','totalnMatrix','totalcMatrix','totalErrorMatrix','cutoff','siteList','n','m')


%%

siteList = [{'dpnmr_larned_east'} {'dpnmr_larned_lwph'} {'dpnmr_larned_west'}...
    {'dpnmr_leque_east'},{'dpnmr_leque_west'}];
load('optimalCutoffTableKansas_var_n_m.mat')
%cutoff = (20:2:800)*10^-3;
 
smoothWindow = 40;
 
larned_eastSmooth = smooth(totalErrorMatrix(1,:),smoothWindow);


smoothErrorMatrix = [g6aSmooth g6bSmooth g5aSmooth g5bSmooth pn1Smooth pn2Smooth g6Smooth g5Smooth];

[minVals, index] = min(smoothErrorMatrix);
%g6aFit = fit(cutoff(15:end)',totalErrorMatrix(15:end,1),'smoothingspline','SmoothingParam',0.99);
cutoff_ms = cutoff*10^3;

minVals = zeros(1, length(siteList));
minIndex = zeros(1, length(siteList));

for j = 1:length(siteList)
   currentErrorCurve = totalErrorMatrix(j,:);
   smoothedError = smooth(currentErrorCurve, smoothWindow);
   
   [minVal, index] = min(smoothedError);
   
   minVals(j) = minVal;
   minIndex(j) = index;
   
   figure(2)
   hold on
   plot(cutoff_ms, currentErrorCurve,'.')
   plot(cutoff_ms, smoothedError)
   plot(cutoff_ms(index),smoothedError(index),'k','MarkerSize',10)
   grid on
   box on
   xlabel('Cutoff (ms')
   ylabel('Mean Average Error (m/day)')
   set(gca,'FontSize',14)
end

legend(siteList)
