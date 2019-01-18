function [K,z,T2dist,T2logbins,k_bootstrap,k_mcmc,k_direct,bestFitMatrix,totalErrorEstimate] = computeProfile(site,n,m,figureson,wDirect)

baseDir = '/Volumes/GoogleDrive/My Drive/USGS Project/USGS Data/';

if strcmp(site,'Site1-WellG5')
    name = 'G5_W1_tr5_20x_16p5_up_F1n2_wRIN_wRFI_Reg50_Va1';
    nmrName = name;
elseif  strcmp(site,'Site1-WellG5above')
    site = 'Site1-WellG5';
    name = 'G5_W1_tr5_20x_16p5_up_F1n2_wRIN_wRFI_Reg50_Va1';
    nmrName = 'G5_W1_tr5_20x_16p5_up_F1n2_wRIN_wRFI_Reg50_Va1_above';
elseif  strcmp(site,'Site1-WellG5below')
    site = 'Site1-WellG5';
    name = 'G5_W1_tr5_20x_16p5_up_F1n2_wRIN_wRFI_Reg50_Va1';
    nmrName = 'G5_W1_tr5_20x_16p5_up_F1n2_wRIN_wRFI_Reg50_Va1_below';
elseif strcmp(site,'Site1-WellG6')
    name = 'G6_W2_tr5_20x_16p75_up_F_wRIN_wRFI_reg50_Va1';
    nmrName = name;
elseif strcmp(site,'Site1-WellG6above')
    site = 'Site1-WellG6';
    name = 'G6_W2_tr5_20x_16p75_up_F_wRIN_wRFI_reg50_Va1';
    nmrName = 'G6_W2_tr5_20x_16p75_up_F_wRIN_wRFI_reg50_Va1_above'
elseif strcmp(site,'Site1-WellG6below')
    site = 'Site1-WellG6';
    name = 'G6_W2_tr5_20x_16p75_up_F_wRIN_wRFI_reg50_Va1';
    nmrName = 'G6_W2_tr5_20x_16p75_up_F_wRIN_wRFI_reg50_Va1_below'
elseif strcmp(site,'Site2-WellPN1')
    name = 'Pl_W1_Tr5_20x_MPp75aLS_F1n2_wRIN_wRFI_Reg50_Va1';
    nmrName = name;
elseif strcmp(site,'Site2-WellPN2')
    name = 'W2_Tr5_20x_MPp75aLS_Reg50_wRIN_wRFI_Va1';
    nmrName = name;
end

in1 = [baseDir site '/' name '/' name '_T2_dist' '.txt']; 
in2 = [baseDir site '/' name '/' name '_T2_bins_log10s' '.txt']; 

T2dist = load(in1); 
T2logbins = load(in2);

[d, K, T2ML, phi, z, SumEch, logK, logT2ML, logPhi, SumEch_3s, SumEch_twm, ...
SumEch_twm_3s] = loadnmrdata2(nmrName); 

logSumEch = log10(SumEch);

bestFitMatrix = zeros(3,3);
totalErrorEstimate = zeros(1,3);


%% Bootstrap

Nboot =  2000; % number of bootstrap samples

% Takes [log10(T2ML), log10(K)] or [log10(T2ML), log10(phi), log10(K)] as a
% single matrix
% [b_boot, n_boot, m_boot] = bootstrap_fun([lt,lp, kk], Nboot);         % m, n can vary
% [b_boot, n_boot, m_boot] = bootstrap_fun([lt, lp, kk], Nboot, n);        % m can vary
%  [b_boot, n_boot, m_boot] = bootstrap_fun_mb([logT2ML, logK], Nboot);    % n can vary
  
%  [b_boot, n_boot, m_boot] = bootstrap_fun([logT2ML, logPhi, logK], Nboot);    % n and m can vary

if isempty(n) && isempty(m)
    [b_boot, n_boot, m_boot] = bootstrap_fun([logT2ML, logPhi, logK], Nboot);
else   
    [b_boot, n_boot, m_boot] = bootstrap_fun([logT2ML, logPhi, logK], Nboot, n, m);   % m, n fixed
end

if figureson ==1 
    bs = log10(b_boot); 
    graph_correlations([bs, n_boot], 2, {'log_{10}(b)', 'n'}, 1, 0)
else
    bs = log10(b_boot); 
end


meanb = mean(b_boot);
sortb = sort(b_boot); 
blo = sortb(50);
bhi = sortb(1950);

median_b = median(b_boot);

% BIG Q: Use mean or median? Maybe should be using median here... 

% Now compute k values as a function of depth
if isempty(m_boot) && isempty(n_boot)
    k_bootstrap = median_b*(phi.^m).*(T2ML).^n;
    
    bestFitMatrix(1,1) = median_b;
    bestFitMatrix(2,1) = m;
    bestFitMatrix(3,1) = n;
    
    totalErrorEstimate(1) = computeError(K, k_bootstrap);
else
    if isempty(m_boot)
        median_n = median(n_boot);
        k_bootstrap = median_b*(phi.^m).*(T2ML).^median_n;
       
        bestFitMatrix(1,1) = median_b;
        bestFitMatrix(2,1) = m;
        bestFitMatrix(3,1) = median_n;
        
        totalErrorEstimate(1) = computeError(K, k_bootstrap);

    else
        median_m = median(m_boot);
        median_n = median(n_boot);
        k_bootstrap = median_b*(phi.^median_m).*(T2ML).^median_n;
        
        bestFitMatrix(1,1) = median_b;
        bestFitMatrix(2,1) = median_m;
        bestFitMatrix(3,1) = median_n;
        
        totalErrorEstimate(1) = computeError(K, k_bootstrap);
    end
end
%% Basic solving for b for fixed n, m
    % Given m and n, we can solve directly for b -> b = log(k) - m*log(phi) -
    % n*log(T2ML). This is the 'direct' method.

if wDirect == 1
    C = @(m, n, lt, lp) m*lp + n*lt; 

    bdir_n2 = logK - C(m, n, logT2ML, logPhi); 

    logb_mean = mean(bdir_n2);
    b_mean = 10.^logb_mean;

    logb_median = median(bdir_n2);
    median_b = 10.^logb_median;

    k_direct = b_mean*(phi.^m).*(T2ML).^2;

    bestFitMatrix(1,2) = median_b;
    bestFitMatrix(2,2) = m;
    bestFitMatrix(3,2) = n;

    totalErrorEstimate(2) = computeError(K, k_direct);
end

%% MCMC for solution to various parameters
% Markov Chain Monte Carlo using Metropolis-Hastings algorithm. Assumes
% Bayes' theorem. Parameters Niter (number of iterations) and stepsize may
% need to be adjusted for good covergence, although with this data set the
% values shown perform reasonably well.
Niter= 1e6; 
stepsize = 0.8; 

%%%%%%%%%%%%%%%%%%%% MCMC hps gives summary parameters (average, CI, etc.) Pick
% either this one or simpler MCMC algorithm to compare with other results.
% NOTE: if the first MCMC algorthm is used, remember when comparing results
% between different methods that the others use fixed values for m and n,
% while this allows both to vary. This will inevitably affect the range in
% b obtained. 

if isempty(m) && isempty(n)

    %%%%%%%%%%%%%%%%%%% MCMC over all parameters (T2B, b, n, m, data error (sigma_mcmc)
    [paramhats, hps, likes, kpreds, accept_rat] = mcmc_nmr_full(K, T2ML, phi, ...
        z, Niter, stepsize, figureson);
    T2B_mcmc  = paramhats(1,:);
    blog_mcmc = paramhats(2,:); 
    n_mcmc = paramhats(3,:);
    m_mcmc = paramhats(4,:);
    sig_mcmc = paramhats(5,:);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Compute b statistics from variable n and m data
    b_mcmc = 10.^blog_mcmc;
    b_mean = mean(b_mcmc);
    b_median = median(b_mcmc);
    
    n_mean = mean(n_mcmc);
    n_median = median(n_mcmc);
    
    m_mean = mean(m_mcmc);
    m_median = median(m_mcmc);
    
    k_mcmc = b_median*(phi.^m_median).*(T2ML).^n_median;
     
    bestFitMatrix(1,3) = b_median;
    bestFitMatrix(2,3) = m_median;
    bestFitMatrix(3,3) = n_median;

else 

    % %%%%%%%%%%%%%% MCMC for b and data error only. n and m is fixed 
%     [paramhats, hps, likes, kpreds, accept_rat] = mcmc_nmr_full(K, T2ML, phi, ...
%         z, Niter, stepsize, figureson);
%     T2B_mcmc  = paramhats(1,:);
%     blog_mcmc = paramhats(2,:); 
%     n_mcmc = paramhats(3,:);
%     m_mcmc = paramhats(4,:);
%     sig_mcmc = paramhats(5,:);
    
    
    [blog_mcmc, sig_mcmc, likes, accept_rat] = mcmc_nmr_bsig (K, T2ML, phi, z, m, n, ...
       Niter, stepsize, figureson); 

    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Compute b statistics from fixed n and m data
    b_mcmc = 10.^blog_mcmc;
    b_mean = mean(b_mcmc);
    b_median = median(b_mcmc);

    k_mcmc = b_median*(phi.^m).*(T2ML).^n;

    bestFitMatrix(1,3) = b_median;
    bestFitMatrix(2,3) = m;
    bestFitMatrix(3,3) = n;

    totalErrorEstimate(3) = computeError(K, k_mcmc);

    params = [blog_mcmc; sig_mcmc];
   %  graph_correlations([bs, n_boot], 2, {'log_{10}(b)', 'n'}, 1, 0)
    graph_correlations(params,2,{'MCMC log_{10}(b)','MCMC \sigma'},1,0)
   
end

% 
% if figureson == 1
%     graph_correlations(paramhats, 2, {'T_B', 'log_{10}(b)', 'n', 'm', '\sigma'}, 0, 0)
%     all_lKpreds = zeros(length(T2ML), length(blog_mcmc));   
% %     for k = 1:length(blog_mcmc)
% %         bbm = [blog_mcmc(k), sig_mcmc(k)];
% %         [~,all_lKpreds(:,k)] = NMRfun2(bbm, K, phi, T2ML, m, n);
% %     end
% %     figure;
% %     hold on
% %     for k = 1:size(all_lKpreds,1)
% %         dpk = all_lKpreds(k,:);
% %         [h,x] = hist(dpk, 40);
% %         t2m = repmat(log10(T2ML(k)), 1, length(x));
% %         scatter(t2m, x, [], h);
% %     end
% %     scatter(log10(T2ML), logK, 'ok', 'filled')
% end

%m_mcmc = 0;


end

