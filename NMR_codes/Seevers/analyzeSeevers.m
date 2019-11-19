% Plot data for Seevers analysis

baseDir = '/Volumes/GoogleDrive/My Drive/Stanford/USGS Project/Field Data/USGS Data/';

load('sonicCoreT2B.mat','sonicCoreT2BData')
T2B_depth = sonicCoreT2BData.Depthm;
T2B_depth = flipud(T2B_depth);

T2B_peak = sonicCoreT2BData.T2Bpeak*10^-3; % Put in s 
T2B_peak = flipud(T2B_peak); 

site = 'Site1-WellG6';
n = 2;
m = 0;
figureson = 1;

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

k_bootstrap = [];
k_mcmc = [];
k_direct = [];

% Going to need to interpolate T2B data and sample at same locations as T2B
% measurements in the lab

% Now compute depths closest to T2B samples
nmrDepths = d(1:end,1);
interpT2B = interp1(T2B_depth, T2B_peak, nmrDepths);

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
seeversT2 = (T2ML.^-1 - interpT2B.^-1).^(-1);
logSeeversT2 = log10(seeversT2);

Nboot = 1000;

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

%%
 %%%%%%%%%%%%%%%%%%% MCMC over all parameters (T2B, b, n, m, data error (sigma_mcmc)
Niter= 1e6; 
stepsize = 0.8; 
 
[paramhats, hps, likes, kpreds, accept_rat] = mcmc_nmr_T2B(K, T2ML, interpT2B, phi, ...
    z, Niter, stepsize, figureson);
blog_mcmc = paramhats(1,:); 
n_mcmc = paramhats(2,:);
m_mcmc = paramhats(3,:);
sig_mcmc = paramhats(4,:);

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

totalErrorEstimate(3) = computeError(K, k_mcmc);

%% Plot data
close all

figure
hold on
plot(T2B_peak,T2B_depth,'*')
plot(interpT2B,z,'+')
set(gca,'YDir','reverse')
legend('T2B','Interp T2B')
grid on
box on


figure
hold on
plot(T2ML,z,'*')
plot(interpT2B/3,z,'+')
set(gca,'YDir','reverse')
legend('T2ML','1/3 T2B')
grid on
box on

T2Seevers = ((T2ML.^(-1) - interpT2B.^(-1)).^(-1));

% Try Fixing T2B
T2Bfix = 500*10^-3;
T2FixedSeevers = ((T2ML.^(-1) - T2Bfix^(-1)).^(-1));

% Try crude optimal T2B
T2Bopt = ones(length(interpT2B),1)*2.2;
T2Bopt(K > 2*10^-4) = 0.520;
T2OptSeevers = ((T2ML.^(-1) - T2Bopt.^(-1)).^(-1));

figure 
hold on
plot(K,T2ML,'*')
%plot(K,T2Seevers,'+')
%plot(K,T2FixedSeevers,'o')
plot(K,T2OptSeevers,'^')
grid 
legend('T2ML','T2Seevers','T2B Fixed')
box on


