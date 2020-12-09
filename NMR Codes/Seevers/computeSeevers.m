x% Compute Seevers Models
%baseDir = '/Volumes/GoogleDrive/My Drive/Stanford/USGS Project/Field Data/USGS Data/';
% baseDir = 'I:\My Drive\Stanford\USGS Project\Field Data\USGS Data\';

%baseDir = '/Volumes/GoogleDrive/My Drive/Stanford/USGS Project/Kansas_Wash_Data/';

[T2dist, T2logbins,nmrName] = loadRawNMRdata(site);

[d, K, T2ML, phi, z, SumEch, logK, logT2ML, logPhi, SumEch_3s, SumEch_twm, ...
SumEch_twm_3s] = loadnmrdata2(nmrName); 

KcutoffLower = 0;
KcutoffHigher = 1*10^(5);

K = K(K>KcutoffLower & K<KcutoffHigher);
phi = phi(K>KcutoffLower & K<KcutoffHigher);
T2ML = T2ML(K>KcutoffLower & K<KcutoffHigher);
SumEch = SumEch(K>KcutoffLower & K<KcutoffHigher);
logK = logK(K>KcutoffLower & K<KcutoffHigher);
logT2ML = logT2ML(K>KcutoffLower & K<KcutoffHigher);
logPhi = logPhi(K>KcutoffLower & K<KcutoffHigher);

logSumEch = log10(SumEch);

bestFitMatrix = zeros(3,3);
totalErrorEstimate = zeros(1,3);

bestFitMatrix(:,2) = NaN;
totalErrorEstimate(2) = NaN;

k_bootstrap = [];
k_mcmc = [];
k_direct = [];

% Going to need to interpolate T2B data and sample at same locations as T2B
% measurements in the lab

% % Now compute depths closest to T2B samples
nmrDepths = d(1:end,1);
%interpT2B = interp1(T2B_depth, T2B_peak, nmrDepths)*10^-3; %Convert from ms to s

%T2B_peak_filt = T2B_peak(T2B_peak > 1800)*10^-3;
%T2Bavg = mean(T2B_peak_filt);
%T2Bavg = 2.2489; %Using a T2B = 2.31s for Maurer data for a GW Temp ~12.5 deg C
T2Bavg = 2.0045;

% Fixing T2B for Kansas/Wash data following Maurer 
%T2Bavg = 2.2;

%% Plot T2B data for reference
% figure
% hold on
% plot(interpT2B, nmrDepths, '*')
% plot(T2B_peak, T2B_depth, '+')
% set(gca,'YDir','reverse')
% legend('Interpolated','Measured')
% grid on
% box on
% xlabel('T2B (s)')
% ylabel('Depth (m)')
% set(gca,'FontSize',14)

%% Implement Bootstrap
%seeversT2 = (T2ML.^-1 - interpT2B.^-1).^(-1);

seeversT2 = (T2ML.^(-1) - T2Bavg.^(-1)).^(-1);
logSeeversT2 = log10(seeversT2);

Nboot = 2000;

if isempty(n) && isempty(m)
    [b_boot, n_boot, m_boot] = bootstrap_fun([logSeeversT2, logPhi, logK], Nboot);
else   
    [b_boot, n_boot, m_boot] = bootstrap_fun([logSeeversT2, logPhi, logK], Nboot, n, m);   % m, n fixed
end

if figureson ==1 
    bs = log10(b_boot); 
%    graph_correlations([bs, n_boot], 2, {'log_{10}(b)', 'n'}, 1, 0)
else
    bs = log10(b_boot); 
end

meanb = mean(b_boot);

median_b = median(b_boot);

% Now compute k values as a function of depth
if isempty(m_boot) && isempty(n_boot)
    k_bootstrap = median_b*(phi.^m).*((T2ML.^(-1)-T2Bavg.^(-1)).^(-1)).^n;

    %k_bootstrap = median_b*(phi.^m).*((T2ML.^(-1)-interpT2B.^(-1)).^(-1)).^n;
    
    bestFitMatrix(1,1) = median_b;
    bestFitMatrix(2,1) = m;
    bestFitMatrix(3,1) = n;
    
    totalErrorEstimate(1) = computeError(K, k_bootstrap);
    totalErrorEstimate(2) = mean(estimateKdiffFactor(K, k_bootstrap, 1));

else
    if isempty(m_boot)
        median_n = median(n_boot);
        k_bootstrap = median_b*(phi.^m).*((T2ML.^(-1)-T2Bavg.^(-1)).^(-1)).^n;
        %k_bootstrap = median_b*(phi.^m).*((T2ML.^(-1)-interpT2B.^(-1)).^(-1)).^n;

        bestFitMatrix(1,1) = median_b;
        bestFitMatrix(2,1) = m;
        bestFitMatrix(3,1) = median_n;
        
        totalErrorEstimate(1) = computeError(K, k_bootstrap);
        totalErrorEstimate(2) = mean(estimateKdiffFactor(K, k_bootstrap, 1));

    else
        median_m = median(m_boot);
        median_n = median(n_boot);
        k_bootstrap = median_b*(phi.^m).*((T2ML.^(-1)-T2Bavg.^(-1)).^(-1)).^n;
        %k_bootstrap = median_b*(phi.^m).*((T2ML.^(-1)-interpT2B.^(-1)).^(-1)).^n;

        
        bestFitMatrix(1,1) = median_b;
        bestFitMatrix(2,1) = median_m;
        bestFitMatrix(3,1) = median_n;
        
        totalErrorEstimate(1) = computeError(K, k_bootstrap);
        totalErrorEstimate(2) = mean(estimateKdiffFactor(K, k_bootstrap, 1));
        
    end
end

%%
%  %%%%%%%%%%%%%%%%%%% MCMC over all parameters (T2B, b, n, m, data error (sigma_mcmc)
% Niter= 1e6; 
% stepsize = 0.8; 
%  
% [paramhats, hps, likes, kpreds, accept_rat] = mcmc_nmr_T2B(K, T2ML, interpT2B, phi, ...
%     z, Niter, stepsize, figureson);
% blog_mcmc = paramhats(1,:); 
% n_mcmc = paramhats(2,:);
% m_mcmc = paramhats(3,:);
% sig_mcmc = paramhats(4,:);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % Compute b statistics from variable n and m data
% b_mcmc = 10.^blog_mcmc;
% b_mean = mean(b_mcmc);
% b_median = median(b_mcmc);
% 
% n_mean = mean(n_mcmc);
% n_median = median(n_mcmc);
% 
% m_mean = mean(m_mcmc);
% m_median = median(m_mcmc);
% 
% k_mcmc = b_median*(phi.^m_median).*(T2ML).^n_median;
% 
% bestFitMatrix(1,3) = b_median;
% bestFitMatrix(2,3) = m_median;
% bestFitMatrix(3,3) = n_median;
% 
% totalErrorEstimate(3) = computeError(K, k_mcmc);




end