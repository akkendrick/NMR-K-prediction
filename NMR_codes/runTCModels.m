% Run Timur-Coates Model Estimates

% Range over pairs of m and n values
close all
clear

siteList = [{'Site1-WellG6'} {'Site1-WellG6above'} {'Site1-WellG6below'}...
    {'Site1-WellG5'} {'Site1-WellG5above'} {'Site1-WellG5below'} {'Site2-WellPN1'} {'Site2-WellPN2'}];

% Matrix Structure... Repeat for each site
%       Bootstrap           
% c
% m
% n

cutoff = 33*10^-3;
m = [0 0 1 1 2 2];
n = [1 2 1 2 1 2];

figureson = 0;
wDirect = 1;
currentRow = -3;

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

        for i = 1:length(siteList)
            baseName = siteList{i}

            [K,z,T2dist,T2logbins,kTC_best,bestFitMatrix,totalErrorEstimate] = computeTCperm(baseName,n(j),m(j),cutoff,figureson);

            currentFitMatrix = [currentFitMatrix bestFitMatrix];
            currentErrorMatrix = [currentErrorMatrix totalErrorEstimate];
        end

        matrixKey(1,j) = n(j);
        matrixKey(2,j) = m(j);
        totalFitMatrix(currentRow:currentRow+2,:) = currentFitMatrix;
        totalErrorMatrix(currentRow,:) = currentErrorMatrix;
    end
end
        