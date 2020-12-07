close all
clear

% Look at width of b distribution

%Wisconsin Data
siteList = [{'Site1-WellG5'}...
   {'Site1-WellG6'} {'Site2-WellPN1'} {'Site2-WellPN2'}];
 
%siteList = [{'Site1-WellG6'}]

% Maurer and Knight
%  siteList = {'dpnmr_larned_east','dpnmr_larned_lwph','dpnmr_larned_west',...
%    'dpnmrA11','dpnmrA12','dpnmrC1S','dpnmrC1SE','dpnmrC1SW',...
%    'dpnmr_leque_east','dpnmr_leque_west'};

m = 4;
n = 2;

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
totalmMatrix = zeros(3,length(m),length(siteList));
totalErrorMatrix = zeros(3,length(m),length(siteList));

tic

 for j = 1:length(m)
    currentRow = currentRow + 4;
    currentFitMatrix = [];
    currentErrorMatrix = [];

    tempb = zeros(3,length(siteList));
    tempm = zeros(3,length(siteList));
    tempError = zeros(3,length(siteList));

   for i = 1:length(siteList)
        siteName = siteList{i};

        if m(j) == 99
            [K,z,k_boot,k_mcmc,k_direct,bestFitMatrix,b_boot,totalErrorEstimate] = computeProfile(siteName,n(j),[],figureson,wDirect,saveData);
        else
            [K,z,k_boot,k_mcmc,k_direct,bestFitMatrix,b_boot,totalErrorEstimate] = computeProfile(siteName,n(j),m(j),figureson,wDirect,saveData);
        end
 
        tempb(:,i) = bestFitMatrix(1,:)';
        tempm(:,i) = bestFitMatrix(2,:)';
        tempError(:,i) = totalErrorEstimate(1,:)';
        
        b_boot_all{i} = b_boot;

    end

    matrixKey(1,j) = n(j);
    matrixKey(2,j) = m(j);

    totalmMatrix(:,j,:) = tempm;
    totalbMatrix(:,j,:) = tempb;
    totalErrorMatrix(:,j,:) = tempError;
end
   
save('SDR_wisc_bboot_n2_m4.mat','siteList','b_boot_all','totalmMatrix','totalbMatrix','totalErrorMatrix')
toc
%%
load('wisc_bboot_n2_m1.mat')

kk = 4;

histogram(b_boot_all{kk});
hold on
grid on
box on

title(siteList{kk})
xlabel('SDR b')
ylabel('Counts')