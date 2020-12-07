
%%
clear

% siteList = [{'Site1-WellG6'} {'Site1-WellG6above'} {'Site1-WellG6below'}...
%     {'Site1-WellG5'} {'Site1-WellG5above'} {'Site1-WellG5below'} {'Site2-WellPN1'} {'Site2-WellPN2'}];

siteList = [{'Site1-WellG6above'} {'Site1-WellG6below'}...
     {'Site1-WellG5above'} {'Site1-WellG5below'} {'Site2-WellPN1'} {'Site2-WellPN2'}];


load 'SDR_allSites_wSep_results_m0_n2.mat'

SDR_Bootstrap = [];
SDR_Direct = [];
SDR_MCMC = [];

% Matrix Structure:
% totalErrorMatrix(model,param combo, site)

for kk = 1:length(siteList)
     [minVal, minbstrpIndex] = min(totalErrorMatrix(1,:,kk));
     SDR_Bootstrap(1,kk) = totalbMatrix(1,minbstrpIndex,kk);
     SDR_Bootstrap(2,kk) = matrixKey(1,minbstrpIndex);
     SDR_Bootstrap(3,kk) = matrixKey(2,minbstrpIndex);
     SDR_Bootstrap(4,kk) = totalErrorMatrix(1,minbstrpIndex,kk);
     
     [minVal, mindirectIndex] = min(totalErrorMatrix(2,:,kk));
     SDR_Direct(1,kk) = totalbMatrix(2,mindirectIndex,kk);
     SDR_Direct(2,kk) = matrixKey(1,mindirectIndex);
     SDR_Direct(3,kk) = matrixKey(2,mindirectIndex);
     SDR_Direct(4,kk) = totalErrorMatrix(2,mindirectIndex,kk);
     
     [minVal, minMCMCIndex] = min(totalErrorMatrix(3,:,kk));
     SDR_MCMC(1,kk) = totalbMatrix(3,minMCMCIndex,kk);
     SDR_MCMC(2,kk) = matrixKey(1,minMCMCIndex);
     SDR_MCMC(3,kk) = matrixKey(2,minMCMCIndex);
     SDR_MCMC(4,kk) = totalErrorMatrix(3,minMCMCIndex,kk);    
end

save('SDR_allSites_wSep_pair_table.mat')

%%
clear

siteList = [{'Site1-WellG6'} {'Site1-WellG6above'} {'Site1-WellG6below'}...
    {'Site1-WellG5'} {'Site1-WellG5above'} {'Site1-WellG5below'} {'Site2-WellPN1'} {'Site2-WellPN2'}];

load 'Seevers_model_pair_results.mat'

Seevers_Bootstrap = [];
Seevers_Direct = [];
Seevers_MCMC = [];

for kk = 1:length(siteList)
     [minVal, minbstrpIndex] = min(totalErrorMatrix(1,:,kk));
     Seevers_Bootstrap(1,kk) = totalbMatrix(1,minbstrpIndex,kk);
     Seevers_Bootstrap(2,kk) = matrixKey(1,minbstrpIndex);
     Seevers_Bootstrap(3,kk) = matrixKey(2,minbstrpIndex);
     Seevers_Bootstrap(4,kk) = totalErrorMatrix(1,minbstrpIndex,kk);
     
     [minVal, mindirectIndex] = min(totalErrorMatrix(2,:,kk));
     Seevers_Direct(1,kk) = totalbMatrix(2,mindirectIndex,kk);
     Seevers_Direct(2,kk) = matrixKey(1,mindirectIndex);
     Seevers_Direct(3,kk) = matrixKey(2,mindirectIndex);
     Seevers_Direct(4,kk) = totalErrorMatrix(2,mindirectIndex,kk);
     
     [minVal, minMCMCIndex] = min(totalErrorMatrix(3,:,kk));
     Seevers_MCMC(1,kk) = totalbMatrix(3,minMCMCIndex,kk);
     Seevers_MCMC(2,kk) = matrixKey(1,minMCMCIndex);
     Seevers_MCMC(3,kk) = matrixKey(2,minMCMCIndex);
     Seevers_MCMC(4,kk) = totalErrorMatrix(3,minMCMCIndex,kk);    
end

save('Seevers_bestFit_table.mat')

%%
clear

siteList = [{'Site1-WellG6'} {'Site1-WellG6above'} {'Site1-WellG6below'}...
    {'Site1-WellG5'} {'Site1-WellG5above'} {'Site1-WellG5below'} {'Site2-WellPN1'} {'Site2-WellPN2'}];

load 'TC_model_results_240.mat'

TC_Bootstrap = [];
TC_Direct = [];
TC_MCMC = [];

for kk = 1:length(siteList)
     [minVal, minbstrpIndex] = min(totalErrorMatrix(1,:,kk));
     TC_Bootstrap(1,kk) = totalbMatrix(1,minbstrpIndex,kk);
     TC_Bootstrap(2,kk) = matrixKey(1,minbstrpIndex);
     TC_Bootstrap(3,kk) = matrixKey(2,minbstrpIndex);
     TC_Bootstrap(4,kk) = totalErrorMatrix(1,minbstrpIndex,kk);
     
     TC_indexQuotient{kk} = totalIndexQuotient{minbstrpIndex}{kk};
     
%      [minVal, mindirectIndex] = min(totalErrorMatrix(2,:,kk));
%      Seevers_Direct(1,kk) = totalbMatrix(2,mindirectIndex,kk);
%      Seevers_Direct(2,kk) = matrixKey(1,mindirectIndex);
%      Seevers_Direct(3,kk) = matrixKey(2,mindirectIndex);
%      Seevers_Direct(4,kk) = totalErrorMatrix(2,mindirectIndex,kk);
%      
%      [minVal, minMCMCIndex] = min(totalErrorMatrix(3,:,kk));
%      Seevers_MCMC(1,kk) = totalbMatrix(3,minMCMCIndex,kk);
%      Seevers_MCMC(2,kk) = matrixKey(1,minMCMCIndex);
%      Seevers_MCMC(3,kk) = matrixKey(2,minMCMCIndex);
%      Seevers_MCMC(4,kk) = totalErrorMatrix(3,minMCMCIndex,kk);    
end

save('TC_bestFit_240_table.mat')

%%
%%
clear

siteList = [{'Site1-WellG6'} {'Site1-WellG6above'} {'Site1-WellG6below'}...
    {'Site1-WellG5'} {'Site1-WellG5above'} {'Site1-WellG5below'} {'Site2-WellPN1'} {'Site2-WellPN2'}];

load 'SDR_model_pair_results.mat'

SDR_Bootstrap = [];
SDR_Direct = [];
SDR_MCMC = [];

% Matrix Structure:
% totalErrorMatrix(model,param combo, site)

for kk = 1:length(siteList)
     [minVal, minbstrpIndex] = min(totalErrorMatrix(1,1:4,kk));
     SDR_Bootstrap(1,kk) = totalbMatrix(1,minbstrpIndex,kk);
     SDR_Bootstrap(2,kk) = matrixKey(1,minbstrpIndex);
     SDR_Bootstrap(3,kk) = matrixKey(2,minbstrpIndex);
     SDR_Bootstrap(4,kk) = totalErrorMatrix(1,minbstrpIndex,kk);
     
     [minVal, mindirectIndex] = min(totalErrorMatrix(2,1:4,kk));
     SDR_Direct(1,kk) = totalbMatrix(2,mindirectIndex,kk);
     SDR_Direct(2,kk) = matrixKey(1,mindirectIndex);
     SDR_Direct(3,kk) = matrixKey(2,mindirectIndex);
     SDR_Direct(4,kk) = totalErrorMatrix(2,mindirectIndex,kk);
     
     [minVal, minMCMCIndex] = min(totalErrorMatrix(3,1:4,kk));
     SDR_MCMC(1,kk) = totalbMatrix(3,minMCMCIndex,kk);
     SDR_MCMC(2,kk) = matrixKey(1,minMCMCIndex);
     SDR_MCMC(3,kk) = matrixKey(2,minMCMCIndex);
     SDR_MCMC(4,kk) = totalErrorMatrix(3,minMCMCIndex,kk);    
end

save('SDR_bestFit_table_for_n1.mat')

%%
clear

siteList = [{'Site1-WellG6'} {'Site1-WellG6above'} {'Site1-WellG6below'}...
    {'Site1-WellG5'} {'Site1-WellG5above'} {'Site1-WellG5below'} {'Site2-WellPN1'} {'Site2-WellPN2'}];

load 'SDR_model_pair_results.mat'

SDR_Bootstrap = [];
SDR_Direct = [];
SDR_MCMC = [];

% Matrix Structure:
% totalErrorMatrix(model,param combo, site)

for kk = 1:length(siteList)
     [minVal, minbstrpIndex] = min(totalErrorMatrix(1,5:8,kk));
     minbstrpIndex = minbstrpIndex + 4;
     
     SDR_Bootstrap(1,kk) = totalbMatrix(1,minbstrpIndex,kk);
     SDR_Bootstrap(2,kk) = matrixKey(1,minbstrpIndex);
     SDR_Bootstrap(3,kk) = matrixKey(2,minbstrpIndex);
     SDR_Bootstrap(4,kk) = totalErrorMatrix(1,minbstrpIndex,kk);
     
     [minVal, mindirectIndex] = min(totalErrorMatrix(2,5:8,kk));
     minbstrpIndex = minbstrpIndex + 4;

     SDR_Direct(1,kk) = totalbMatrix(2,mindirectIndex,kk);
     SDR_Direct(2,kk) = matrixKey(1,mindirectIndex);
     SDR_Direct(3,kk) = matrixKey(2,mindirectIndex);
     SDR_Direct(4,kk) = totalErrorMatrix(2,mindirectIndex,kk);
     
     [minVal, minMCMCIndex] = min(totalErrorMatrix(3,5:8,kk));
     minbstrpIndex = minbstrpIndex + 4;

     SDR_MCMC(1,kk) = totalbMatrix(3,minMCMCIndex,kk);
     SDR_MCMC(2,kk) = matrixKey(1,minMCMCIndex);
     SDR_MCMC(3,kk) = matrixKey(2,minMCMCIndex);
     SDR_MCMC(4,kk) = totalErrorMatrix(3,minMCMCIndex,kk);    
end

save('SDR_bestFit_table_for_n2.mat')
