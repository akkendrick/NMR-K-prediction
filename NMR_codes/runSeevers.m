% Run Seevers Models

% Range over pairs of m and n values
close all
clear

siteList = [{'Site1-WellG6'} {'Site1-WellG6above'} {'Site1-WellG6below'}...
    {'Site1-WellG5'} {'Site1-WellG5above'} {'Site1-WellG5below'} {'Site2-WellPN1'} {'Site2-WellPN2'}];

m = [0 1 2 0 1 2];
n = [1 1 1 2 2 2];

figureson = 1;
wDirect = 0;
currentRow = -3;

% % for k = 1:length(siteList)
% % 
% %     computeSeevers(siteList{k},n,m,figureson,0)
% %     title(siteList{k})
% % 
% % end
if isempty(m) && isempty(n)
    currentFitMatrix = [];
    currentErrorMatrix = [];
    for i = 1:length(siteList)
        siteName = siteList{i}

        [K,z,T2dist,T2logbins,k_boot,k_mcmc,k_direct,bestFitMatrix,totalErrorEstimate] = computeSeevers(siteName,[],[],figureson,wDirect);

        currentFitMatrix = [currentFitMatrix bestFitMatrix];
        currentErrorMatrix = [currentErrorMatrix totalErrorEstimate];
    end
else
    for j = 1:length(m)
        currentRow = currentRow + 4;
        currentFitMatrix = [];
        currentErrorMatrix = [];

        for i = 1:length(siteList)
            siteName = siteList{i}

            [K,z,T2dist,T2logbins,k_boot,k_mcmc,k_direct,bestFitMatrix,totalErrorEstimate] = computeSeevers(siteName,n(j),m(j),figureson,wDirect);

            currentFitMatrix = [currentFitMatrix bestFitMatrix];
            currentErrorMatrix = [currentErrorMatrix totalErrorEstimate];
        end

        matrixKey(1,j) = n(j);
        matrixKey(2,j) = m(j);
        totalFitMatrix(currentRow:currentRow+2,:) = currentFitMatrix;
        totalErrorMatrix(currentRow,:) = currentErrorMatrix;
    end
end
