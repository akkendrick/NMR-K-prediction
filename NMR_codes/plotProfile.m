% Plot profile and compare different estimates
clear
close all

site = 'Site1-WellG6below'

%%
n = 2;
m = 0;
figureson = 1;
wDirect = 1;

[K,z,T2dist,T2logbins,k_boot,k_mcmc,k_direct,bestFitMatrix,totalErrorEstimate] = computeProfile(site,n,m,figureson,wDirect);

k_estimates = [k_boot k_direct k_mcmc];
k_names = [{'DPP'} {'Bootstrap n=2 m=4'} {'Direct n=2 m=4'} {'MCMC n=2 m=4'}];
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
%plotKwithDepth(K,z,T2dist,T2logbins,k_estimates,k_names,k_sym)

%%

