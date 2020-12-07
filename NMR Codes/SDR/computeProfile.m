function [K,z,k_bootstrap,k_mcmc,k_direct,bestFitMatrix,b_boot,totalErrorEstimate] = computeProfile(site,n,m,figureson,wDirect,saveData)

baseDir = '/Volumes/GoogleDrive/My Drive/USGS Project/USGS Data/';
%baseDir = 'I:\My Drive\USGS Project\USGS Data\';

if strcmp(site,'Site1-WellG5')
    name = 'G5_W1_tr5_20x_16p5_up_F1n2_wRIN_wRFI_Reg50_Va1'; 
    nmrName = name;
    saveName = 'G5';
   
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
    saveName = 'G6';
    
elseif strcmp(site,'Site1-WellG6above')
    site = 'Site1-WellG6';
    name = 'G6_W2_tr5_20x_16p75_up_F_wRIN_wRFI_reg50_Va1';
    nmrName = 'G6_W2_tr5_20x_16p75_up_F_wRIN_wRFI_reg50_Va1_above';
    
elseif strcmp(site,'Site1-WellG6below')
    site = 'Site1-WellG6';
    name = 'G6_W2_tr5_20x_16p75_up_F_wRIN_wRFI_reg50_Va1';
    nmrName = 'G6_W2_tr5_20x_16p75_up_F_wRIN_wRFI_reg50_Va1_below';
elseif strcmp(site,'Site2-WellPN1')
    name = 'Pl_W1_Tr5_20x_MPp75aLS_F1n2_wRIN_wRFI_Reg50_Va1';
    nmrName = name;
    saveName = 'PN1';
    
elseif strcmp(site,'Site2-WellPN2')
    name = 'W2_Tr5_20x_MPp75aLS_Reg50_wRIN_wRFI_Va1';
    nmrName = name;
    saveName = 'PN2';
    
else
    name = site;
    nmrName = site;
    saveName = site;
end

%if figureson == 1
%     in1 = [baseDir site '/' name '/' name '_T2_dist' '.txt']; 
%     in2 = [baseDir site '/' name '/' name '_T2_bins_log10s' '.txt']; 
% 
%     T2dist = load(in1); 
%     T2logbins = load(in2);
% %end

[d, K, T2ML, phi, z, SumEch, logK, logT2ML, logPhi, SumEch_3s, SumEch_twm, ...
SumEch_twm_3s] = loadnmrdata2(nmrName); 

logSumEch = log10(SumEch);

bestFitMatrix = zeros(3,3);
totalErrorEstimate = zeros(1,3);

k_bootstrap = [];
k_mcmc = [];
k_direct = [];
%% Bootstrap

% Provide bootstrap an evenly weighted distribution points
% Idea:

%edges = [10^-6,5*10^-6,10^-5,5*10^-5,10^-4,5*10^-4,10^-3];
%edges = [10^-6,10^-4,10^-3];


% edges = [10^-6,5*10^-5,10^-3];
% discretizedK = discretize(K,edges);
% 
% nPts = length(K);
% %weights(discretizedK == 1) =  1-sum(discretizedK == 1)/nPts;
% % weights(discretizedK == 1) =  1;
% % 
% % weights(discretizedK == 2) =  1-sum(discretizedK == 2)/[[nPts;
% 
% [targetPoints, index] = min([sum(discretizedK == 1), sum(discretizedK == 2)]);
% 
% sampledK_2 = datasample(K(discretizedK == 2),targetPoints,'Replace',false); 
% sampledK_1 = datasample(K(discretizedK == 1),targetPoints,'Replace',false);
% 
% sampledK = [sampledK_2; sampledK_1];
% 
% discretizedSample = discretize(sampledK,edges);
% hist(discretizedSample)
% 
% % Now sampledK has the small K and a random selection of high K samples
% for kk = 1:length(sampledK)
%     % Take the min just in case we have multiple K values that are the same
%     filtIndex = min(find(K == sampledK(kk)));
%     logT2ML_filt(kk) = logT2ML(filtIndex);
%     logPhi_filt(kk) = logPhi(filtIndex);
%     logK_filt(kk) = logK(filtIndex);
% end
% 
% logT2ML_filt = logT2ML_filt';
% logPhi_filt = logPhi_filt';
% logK_filt = logK_filt';



logT2ML_filt = logT2ML;
logPhi_filt = logPhi;
logK_filt = logK;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Edges for Wisconsin
% Only look at "small values of K"
%lowestKEdge = 6*10^-5;
%midKEdge = 10^-4; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Edges for Maurer and Knight
% Only look at "small values of K"
% lowestKEdge = 10^-5;
% midKEdge = 5*10^-4; 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %filtInd = find(K < lowestKEdge);
% %filtInd = find(K < midKEdge & K > lowestKEdge);
% filtInd = find(K > midKEdge);
% 
% logT2ML_filt = logT2ML(filtInd);
% logPhi_filt = logPhi(filtInd);
% logK_filt = logK(filtInd);

%%
% Use discretize to group data into bins
% Assign each bin a weight based on the relative freq compared to the total
% Use datasample with weights to try and obtain an even distribution, plot
% hist for mult realizations, check results


Nboot =  2000; % number of bootstrap samples

% Takes [log10(T2ML), log10(K)] or [log10(T2ML), log10(phi), log10(K)] as a
% single matrix
% [b_boot, n_boot, m_boot] = bootstrap_fun([lt,lp, kk], Nboot);         % m, n can vary
% [b_boot, n_boot, m_boot] = bootstrap_fun([lt, lp, kk], Nboot, n);        % m can vary
%  [b_boot, n_boot, m_boot] = bootstrap_fun_mb([logT2ML, logK], Nboot);    % n can vary
  
%  [b_boot, n_boot, m_boot] = bootstrap_fun([logT2ML, logPhi, logK], Nboot);    % n and m can vary

% if isempty(n) && isempty(m)
%     [b_boot, n_boot, m_boot] = bootstrap_fun([logT2ML, logPhi, logK], Nboot);
% elseif ~isempty(n) && isempty(m)
%     [b_boot, n_boot, m_boot] = bootstrap_fun([logT2ML, logPhi, logK], Nboot, n);
% else
%     [b_boot, n_boot, m_boot] = bootstrap_fun([logT2ML, logPhi, logK], Nboot, n, m);   % m, n fixed
% end

if isempty(n) && isempty(m)
    [b_boot, n_boot, m_boot] = bootstrap_fun([logT2ML_filt, logPhi_filt, logK_filt], Nboot);
elseif ~isempty(n) && isempty(m)
    [b_boot, n_boot, m_boot] = bootstrap_fun([logT2ML_filt, logPhi_filt, logK_filt], Nboot, n);
else
    [b_boot, n_boot, m_boot] = bootstrap_fun([logT2ML_filt, logPhi_filt, logK_filt], Nboot, n, m);   % m, n fixed
end

if figureson ==1 
    bs = log10(b_boot); 
    graph_correlations([bs, n_boot], 2, {'log_{10}(b)', 'n'}, 1, 0)
else
    bs = log10(b_boot); 
end


meanb = mean(b_boot);
sortb = sort(b_boot); 
% blo = sortb(50);
% bhi = sortb(1950);

median_b = median(b_boot)
std_b = std(b_boot);

% BIG Q: Use mean or median? Maybe should be using median here... 

% Now compute k values as a function of depth
if isempty(m_boot) && isempty(n_boot)
    k_bootstrap = median_b*(phi.^m).*(T2ML).^n;
    
    bestFitMatrix(1,1) = median_b;
    bestFitMatrix(2,1) = m;
    bestFitMatrix(3,1) = n;
    
    totalErrorEstimate(1) = computeError(K, k_bootstrap);
else
    if isempty(n_boot) && ~isempty(m_boot)
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
        totalErrorEstimate(2) = mean(estimateKdiffFactor(K,k_bootstrap,1))
        totalErrorEstimate(3) = std_b;
    end
end

%graph_correlations([b_boot, m_boot], 2, {'log_{10}(b)', 'm'}, 1, 0)

if saveData ==1 
    save(strcat(saveName,'_bootstrap_n_m_var.mat'),'bs','n_boot','m_boot')
end
%% Basic solving for b for fixed n, m
    % Given m and n, we can solve directly for b -> b = log(k) - m*log(phi) -
    % n*log(T2ML). This is the 'direct' method.

if wDirect == 1
    mDirect = bestFitMatrix(2,1);
    nDirect = bestFitMatrix(3,1);
    
    C = @(m, n, lt, lp) m*lp + n*lt; 

    bdir_n2 = logK - C(mDirect, nDirect, logT2ML, logPhi); 

    logb_mean = mean(bdir_n2);
    b_mean = 10.^logb_mean;
    bs_basic = 10.^(bdir_n2);


    logb_median = median(bdir_n2);
    median_b = 10.^logb_median;

    k_direct = b_mean*(phi.^mDirect).*(T2ML).^nDirect;

    bestFitMatrix(1,2) = median_b;
    bestFitMatrix(2,2) = mDirect;
    bestFitMatrix(3,2) = nDirect;

    totalErrorEstimate(2) = computeError(K, k_direct);
end

if saveData == 1
    save(strcat(saveName,'_basic_solving_n_m_var.mat'),'bs_basic','nDirect','mDirect')
end

% %% MCMC for solution to various parameters
% % Markov Chain Monte Carlo using Metropolis-Hastings algorithm. Assumes
% % Bayes' theorem. Parameters Niter (number of iterations) and stepsize may
% % need to be adjusted for good covergence, although with this data set the
% % values shown perform reasonably well.
% Niter= 1e6; 
% stepsize = 0.8; 
% 
% %%%%%%%%%%%%%%%%%%%% MCMC hps gives summary parameters (average, CI, etc.) Pick
% % either this one or simpler MCMC algorithm to compare with other results.
% % NOTE: if the first MCMC algorthm is used, remember when comparing results
% % between different methods that the others use fixed values for m and n,
% % while this allows both to vary. This will inevitably affect the range in
% % b obtained. 
% 
% if isempty(m) || isempty(n)
% 
%     %%%%%%%%%%%%%%%%%%% MCMC over all parameters (T2B, b, n, m, data error (sigma_mcmc)
%     [paramhats, hps, likes, kpreds, accept_rat] = mcmc_nmr_full(K, T2ML, phi, ...
%         z, Niter, stepsize, figureson);
%     T2B_mcmc  = paramhats(1,:);
%     blog_mcmc = paramhats(2,:); 
%     n_mcmc = paramhats(3,:);
%     m_mcmc = paramhats(4,:);
%     sig_mcmc = paramhats(5,:);
% 
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     
%     % Compute b statistics from variable n and m data
%     b_mcmc = 10.^blog_mcmc;
%     b_mean = mean(b_mcmc);
%     b_median = median(b_mcmc);
%     
%     n_mean = mean(n_mcmc);
%     n_median = median(n_mcmc);
%     
%     m_mean = mean(m_mcmc);
%     m_median = median(m_mcmc);
%     
%     k_mcmc = b_median*(phi.^m_median).*(T2ML).^n_median;
%      
%     bestFitMatrix(1,3) = b_median;
%     bestFitMatrix(2,3) = m_median;
%     bestFitMatrix(3,3) = n_median;
% 
%     totalErrorEstimate(3) = computeError(K, k_mcmc);
% 
% else 
% 
%     %%%%%%%%%%%%%% MCMC for b and data error only. n and m is fixed 
% %     [paramhats, hps, likes, kpreds, accept_rat] = mcmc_nmr_full(K, T2ML, phi, ...
% %         z, Niter, stepsize, figureson);
% %     T2B_mcmc  = paramhats(1,:);
% %     blog_mcmc = paramhats(2,:); 
% %     n_mcmc = paramhats(3,:);
% %     m_mcmc = paramhats(4,:);
% %     sig_mcmc = paramhats(5,:);
%     
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%     
%     
%     [blog_mcmc, sig_mcmc, likes, accept_rat] = mcmc_nmr_bsig (K, T2ML, phi, z, m, n, ...
%        Niter, stepsize, figureson); 
% 
%     
%  %   Compute b statistics from fixed n and m data
%     b_mcmc = 10.^blog_mcmc;
%     b_mean = mean(b_mcmc);
%     b_median = median(b_mcmc);
%     
%     k_mcmc = b_median*(phi.^m).*(T2ML).^n;
%     
%     bestFitMatrix(1,3) = b_median;
%     bestFitMatrix(2,3) = m;
%     bestFitMatrix(3,3) = n;
%     
%     totalErrorEstimate(3) = computeError(K, k_mcmc);
% 
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     
%     params = [blog_mcmc; sig_mcmc];
%    %  graph_correlations([bs, n_boot], 2, {'log_{10}(b)', 'n'}, 1, 0)
%     graph_correlations(params,2,{'MCMC log_{10}(b)','MCMC \sigma'},1,0)
%    
%     
%     
% end
% 
% if saveData == 1
%     save(strcat(saveName,'_MCMC_n_m_var.mat'),'blog_mcmc','n_mcmc','m_mcmc')
% end
% 
% % 
% % if figureson == 1
% %     graph_correlations(paramhats, 2, {'T_B', 'log_{10}(b)', 'n', 'm', '\sigma'}, 0, 0)
% %     all_lKpreds = zeros(length(T2ML), length(blog_mcmc));   
% % %     for k = 1:length(blog_mcmc)
% % %         bbm = [blog_mcmc(k), sig_mcmc(k)];
% % %         [~,all_lKpreds(:,k)] = NMRfun2(bbm, K, phi, T2ML, m, n);
% % %     end
% % %     figure;
% % %     hold on
% % %     for k = 1:size(all_lKpreds,1)
% % %         dpk = all_lKpreds(k,:);
% % %         [h,x] = hist(dpk, 40);
% % %         t2m = repmat(log10(T2ML(k)), 1, length(x));
% % %         scatter(t2m, x, [], h);
% % %     end
% % %     scatter(log10(T2ML), logK, 'ok', 'filled')
% % end
% 
% %m_mcmc = 0;


end

