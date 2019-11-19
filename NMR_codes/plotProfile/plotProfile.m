% Plot profile and compare different estimates
clear
close all

site = 'Site1-WellG6'

%%
n = 2;
m = 0;
figureson = 1;
wDirect = 1;

[K,z,T2dist,T2logbins,k_boot,k_mcmc,k_direct,bestFitMatrix,totalErrorEstimate] = computeProfile(site,n,m,figureson,wDirect);

k_estimates = [k_boot k_direct k_mcmc];
k_names = [{'DPP'} {'Bootstrap'} {'Direct'} {'MCMC'}];
k_sym = [{'+'} {'+'} {'+'}];

plotKwithDepth(K,z,T2dist,T2logbins,k_estimates,k_names,k_sym)

methodNames = [{'Bootstrap'} {'Direct'} {'MCMC'}]

bestFitMatrix
totalErrorEstimate

%%
% n = 2;
% m = 1;
% [K,z,T2dist,T2logbins,k_boot,k_mcmc,~,~] = computeProfile(site,n,m);
% 
% k_estimates = [k_boot k_mcmc];
% k_names = [{'DPP'} {'Bootstrap n=2 m=1'} {'MCMC n=2 m=1'}];
% k_sym = [{'+'} {'+'}];
%  
% 
% [K,z,T2dist,T2logbins,k_boot,k_mcmc,~,~] = computeProfile(site,[],[]);
% k_estimates = [k_estimates k_boot k_mcmc];
% k_names = [k_names {'Best Bootstrap'} {'Best MCMC'}];
% k_sym = [k_sym {'o'} {'o'}];
%%
clear
close all

site = 'Site1-WellG5above'

tic

figureson = 0;
wDirect = 1;

n = 2;
m = 0;

[K,z,T2dist,T2logbins,k_boot,k_mcmc,k_direct,bestFitMatrix,totalErrorEstimate] = computeProfile(site,n,m,figureson,wDirect);

k_estimates = [k_boot k_mcmc];
k_names = [{'1:1'} {'Bootstrap n2 m0'} {'MCMC n2 m0'}];
k_sym = [{'+'} {'+'}];

n = 2;
m = 1;

[K,z,T2dist,T2logbins,k_boot,k_mcmc,k_direct,bestFitMatrix,totalErrorEstimate] = computeProfile(site,n,m,figureson,wDirect);

k_estimates = [k_estimates k_boot k_mcmc];
k_names = [k_names {'Bootstrap n2 m1'} {'MCMC n2 m1'}];
k_sym = [k_sym {'*'} {'*'}];

n = 2;
m = 4;

[K,z,T2dist,T2logbins,k_boot,k_mcmc,k_direct,bestFitMatrix,totalErrorEstimate] = computeProfile(site,n,m,figureson,wDirect);

k_estimates = [k_estimates k_boot k_mcmc];
k_names = [k_names {'Bootstrap n2 m4'} {'MCMC n2 m4'}];
k_sym = [k_sym {'o'} {'o'}];


plotKwithDepth(K,z,T2dist,T2logbins,k_estimates,k_names,k_sym)
plotKestKdpp(K,k_estimates,k_names,k_sym)

title(site)

toc

%%

clear
close all

site = 'Site1-WellG6below'

tic

figureson = 0;
wDirect = 1;

n = 2;
m = 1;

[K,z,T2dist,T2logbins,k_boot,k_mcmc,k_direct,bestFitMatrix,totalErrorEstimate] = computeProfile(site,n,m,figureson,wDirect);

k_estimates = [k_boot];
k_names = [{'1:1'} {'Bootstrap n2 m1'} ];
k_sym = [{'+'}];

n = 1;
m = 1;

[K,z,T2dist,T2logbins,k_boot,k_mcmc,k_direct,bestFitMatrix,totalErrorEstimate] = computeProfile(site,n,m,figureson,wDirect);

k_estimates = [k_estimates k_boot];
k_names = [k_names {'Bootstrap n1 m1'}];
k_sym = [k_sym {'*'}];

n = 0;
m = 1;

[K,z,T2dist,T2logbins,k_boot,k_mcmc,k_direct,bestFitMatrix,totalErrorEstimate] = computeProfile(site,n,m,figureson,wDirect);

k_estimates = [k_estimates k_boot];
k_names = [k_names {'Bootstrap n0 m1'}];
k_sym = [k_sym {'o'}];


plotKwithDepth(K,z,T2dist,T2logbins,k_estimates,k_names,k_sym)
plotKestKdpp(K,k_estimates,k_names,k_sym)

title(site)

toc

