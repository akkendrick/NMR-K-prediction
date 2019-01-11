% Plot measured and estimated k with depth

clear

%name = 'G5_W1_tr5_20x_16p5_up_F1n2_wRIN_wRFI_Reg50_Va1';
%name = 'G6_W2_tr5_20x_16p75_up_F_wRIN_wRFI_reg50_Va1';
%name = 'Pl_W1_Tr5_20x_MPp75aLS_F1n2_wRIN_wRFI_Reg50_Va1';
name = 'W2_Tr5_20x_MPp75aLS_Reg50_wRIN_wRFI_Va1'


[d, K, T2ML, phi, z, SumEch, logK, logT2ML, logPhi, SumEch_3s, SumEch_twm, ...
SumEch_twm_3s] = loadnmrdata2(name); 

logSumEch = log10(SumEch); 

% Estimate SDR parameters from bootstrap
%% Bootstrap
% Randomly sample data with replacement, solve the subset for the
% best-fitting parameter values, and repeat many times. 

% m assumed 0; 
n = 2; 
m = 4; 
Nboot =  2000; % number of bootstrap samples

% Takes [log10(T2ML), log10(K)] or [log10(T2ML), log10(phi), log10(K)] as a
% single matrix
% [b_boot, n_boot, m_boot] = bootstrap_fun([lt,lp, kk], Nboot);         % m, n can vary
% [b_boot, n_boot, m_boot] = bootstrap_fun([lt, lp, kk], Nboot, n);        % m can vary
% [b_boot, n_boot, m_boot] = bootstrap_fun_mb([logT2ML, logK], Nboot);    % n can vary
 [b_boot, n_boot, m_boot] = bootstrap_fun([logT2ML, logPhi, logK], Nboot, n, m);   % m, n fixed

bs = log10(b_boot); 
graph_correlations([bs, n_boot], 2, {'log_{10}(b)', 'n'}, 1, 0)

meanb = mean(b_boot)
sortb = sort(b_boot); 
blo = sortb(50)
bhi = sortb(1950)

% Now compute k values as a function of depth

k_bootstrap = meanb*(phi.^m).*(T2ML).^n;

%plotKwithz2(K, z, k_bootstrap)
plotKwithDepth(K,z,k_bootstrap)

%% Basic solving for b for fixed n, m
    % Given m and n, we can solve directly for b -> b = log(k) - m*log(phi) -
    % n*log(T2ML). This is the 'direct' method.

C = @(m, n, lt, lp) m*lp + n*lt; 
n = 2;
m = 4; 
bdir_n1 = logK - C(m, n, logT2ML, logPhi); 

n = 2; 
bdir_n2 = logK - C(m, n, logT2ML, logPhi); 

logb_mean = mean(bdir_n2);
b_mean = 10.^logb_mean;

logb_median = median(bdir_n2);
b_median = 10.^logb_median;

k_direct = b_mean*(phi.^m).*(T2ML).^2;

plotKwithDepth(K,z,k_direct)

%%
