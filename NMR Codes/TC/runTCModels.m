% Run Timur-Coates Model Estimates

% Range over pairs of m and n values
close all
clear

%siteList = [{'Site1-WellG6'} {'Site1-WellG6above'} {'Site1-WellG6below'}...
%    {'Site1-WellG5'} {'Site1-WellG5above'} {'Site1-WellG5below'} {'Site2-WellPN1'} {'Site2-WellPN2'}];

siteList = [{'Site1-WellG6'} {'Site1-WellG5'} {'Site2-WellPN1'} {'Site2-WellPN2'}];

% Matrix Structure... Repeat for each site
%       Bootstrap           
% c
% m
% n

%cutoff = 33*10^-3;
%m = [0 1 2 4 0 1 2 4] ;
%n = [1 1 1 1 2 2 2 2];

cutoff = 42.1*10^-3;
m = [1];
n = [2];

figureson = 0;
wDirect = 1;
currentRow = -3;

tic

if isempty(m) && isempty(n)
    currentFitMatrix = [];
    currentErrorMatrix = [];
    for i = 1:length(siteList)
        baseName = siteList{i}

        [K,z,T2dist,T2logbins,kTC_best,bestFitMatrix,totalErrorEstimate] = computeTCperm(baseName,[],[],cutoff,figureson);

        currentFitMatrix = [currentFitMatrix bestFitMatrix];
        currentErrorMatrix = [currentErrorMatrix totalErrorEstimate];
    end
else
    for j = 1:length(m)
        currentRow = currentRow + 4;
        currentFitMatrix = [];
        currentErrorMatrix = [];

        tempb = zeros(3,length(siteList));
        tempError = zeros(3,length(siteList));
        tempQuotient = {};

        for i = 1:length(siteList)
            baseName = siteList{i};

            [K,z,T2dist,T2logbins,kTC_best,bestFitMatrix,totalErrorEstimate,indexQuotient] = computeTCperm(baseName,n(j),m(j),cutoff,figureson);
             
            plotKestKdpp(K,kTC_best,[{'DPP'} {'TC'}],[{'+'} {'+'}])
            
            tempb(:,i) = bestFitMatrix(1,:)';
            tempError(:,i) = totalErrorEstimate(1,:)';
            tempQuotient{i} = indexQuotient;
        end

        matrixKey(1,j) = n(j);
        matrixKey(2,j) = m(j);
        
        totalbMatrix(:,j,:) = tempb;
        totalErrorMatrix(:,j,:) = tempError;
        totalIndexQuotient{:,j} = tempQuotient;
    end
end

toc
save('TC_results.mat','m','n','matrixKey','totalbMatrix','totalErrorMatrix','totalIndexQuotient')

        