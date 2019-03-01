% Range over pairs of m and n values
close all
clear

siteList = [{'Site1-WellG6'} {'Site1-WellG6above'} {'Site1-WellG6below'}...
    {'Site1-WellG5'} {'Site1-WellG5above'} {'Site1-WellG5below'} {'Site2-WellPN1'} {'Site2-WellPN2'}];

%siteList =  [{'Site1-WellG6'}] %{'Site1-WellG6above'}]

% Matrix Structure... Repeat for each site
%       Bootstrap     Direct      MCMC
% b
% m
% n
%             ROW OF ZEROES

% m = [0 0 2 2 4 4];
% n = [1 2 1 2 1 2];

m = [0 1 2 4 0 1 2 4];
n = [1 1 1 1 2 2 2 2];

% m = [0 1 4]
% n = [0 0 0]

% m = [];
% n = [];

figureson = 0;
wDirect = 1;
currentRow = -3;
saveData = 0;

mcmcMatrix = [];
bootMatrix = [];
directMatrix = [];

% Matrix is organized as follows: 
% (statModel, m and n combo, site)

% For each site, rows correspond to a different statistical model
% Where the stat models are organized as:
%   bootstrap
%   direct
%   mcmc

% Each column is then for a different pair of m and n values

totalbMatrix = zeros(3,length(m),length(siteList));
totalErrorMatrix = zeros(3,length(m),length(siteList));

tic

if isempty(m) && isempty(n)
    currentFitMatrix = [];
    currentErrorMatrix = [];
    for i = 1:length(siteList)
        siteName = siteList{i}

        [K,z,T2dist,T2logbins,k_boot,k_mcmc,k_direct,bestFitMatrix,totalErrorEstimate] = computeProfile(siteName,[],[],figureson,wDirect,saveData);

        currentFitMatrix = [currentFitMatrix bestFitMatrix];
        currentErrorMatrix = [currentErrorMatrix totalErrorEstimate];
    end
    
elseif isempty(m) && ~isempty(n)
    for j = 1:length(n)
        currentRow = currentRow + 4;
        currentFitMatrix = [];
        currentErrorMatrix = [];
        
        tempb = zeros(3,length(siteList));
        tempError = zeros(3,length(siteList));
        
        for i = 1:length(siteList)
            siteName = siteList{i};
            [K,z,T2dist,T2logbins,k_boot,k_mcmc,k_direct,bestFitMatrix,totalErrorEstimate] = computeProfile(siteName,n(j),[],figureson,wDirect,saveData);
            
            tempb(:,i) = bestFitMatrix(1,:)';
            tempError(:,i) = totalErrorEstimate(1,:)';
            
        end

        matrixKey(1,j) = n(j);
        matrixKey(2,j) = m(j);
        
        totalbMatrix(:,j,:) = tempb;
        totalErrorMatrix(:,j,:) = tempError;
        
    end
else
    for j = 1:length(m)
        currentRow = currentRow + 4;
        currentFitMatrix = [];
        currentErrorMatrix = [];
        
        tempb = zeros(3,length(siteList));
        tempError = zeros(3,length(siteList));
        
        for i = 1:length(siteList)
            siteName = siteList{i};
            [K,z,T2dist,T2logbins,k_boot,k_mcmc,k_direct,bestFitMatrix,totalErrorEstimate] = computeProfile(siteName,n(j),m(j),figureson,wDirect,saveData);
            
            tempb(:,i) = bestFitMatrix(1,:)';
            tempError(:,i) = totalErrorEstimate(1,:)';
            
        end

        matrixKey(1,j) = n(j);
        matrixKey(2,j) = m(j);
        
        totalbMatrix(:,j,:) = tempb;
        totalErrorMatrix(:,j,:) = tempError;
        
    end
end
        
toc            

save('SDR_model_pair_results.mat','m','n','matrixKey','totalbMatrix','totalErrorMatrix')
            
            