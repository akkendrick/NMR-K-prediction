% Slow vs Fast Diffusion Comparision

clear
%close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sites = [{'Site1-WellG6'} {'Site1-WellG5'} {'Site2-WellPN1'} {'Site2-WellPN2'}];
figureNames = {'Adams G6','Adams G5','Plainfield PN1','Plainfield PN2'}

load SDR_bestFit_1014_m0_n2.mat

sites = siteList;
SDR_b_fast = squeeze(totalbMatrix(1,1,:))';
SDR_n_fast = n;
SDR_n_fast = repmat(n,length(siteList),1);
SDR_m_fast = squeeze(totalmMatrix(1,1,:))';

load Seevers_bestFit_1014_m0_n2.mat
sites = siteList;
Seevers_b_fast = squeeze(totalbMatrix(1,1,:))';
Seevers_n_fast = n;
Seevers_n_fast = repmat(n,length(siteList),1);
Seevers_m_fast = squeeze(totalmMatrix(1,1,:))';


load SDR_bestFit_1014_m0_n1.mat
sites = siteList;
SDR_b_slow = squeeze(totalbMatrix(1,1,:))';
SDR_n_slow = n;
SDR_n_slow = repmat(n,length(siteList),1);
SDR_m_slow = squeeze(totalmMatrix(1,1,:))';

load Seevers_bestFit_1014_m0_n1.mat
sites = siteList;
Seevers_b_slow = squeeze(totalbMatrix(1,1,:))';
Seevers_n_slow = n;
Seevers_n_slow = repmat(n,length(siteList),1);
Seevers_m_slow = squeeze(totalmMatrix(1,1,:))';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sites = {'dpnmr_larned_east','dpnmr_larned_lwph','dpnmr_larned_west',...
%   'dpnmrA11','dpnmrA12','dpnmrC1S','dpnmrC1SE','dpnmrC1SW',...
%   'dpnmr_leque_east','dpnmr_leque_west'};
% 
% load SDR_maurer_bestFit_1015_m0_n2.mat
%     
% sites = siteList;
% SDR_b_fast = squeeze(totalbMatrix(1,1,:))';
% SDR_n_fast = n;
% SDR_n_fast = repmat(n,length(siteList),1);
% SDR_m_fast = squeeze(totalmMatrix(1,1,:))';
% 
% load Seevers_maurer_bestFit_1015_m0_n2.mat
% sites = siteList;
% Seevers_b_fast = squeeze(totalbMatrix(1,1,:))';
% Seevers_n_fast = n;
% Seevers_n_fast = repmat(n,length(siteList),1);
% Seevers_m_fast = squeeze(totalmMatrix(1,1,:))';
% 
% 
% load SDR_maurer_bestFit_1015_m0_n1.mat
% sites = siteList;
% SDR_b_slow = squeeze(totalbMatrix(1,1,:))';
% SDR_n_slow = n;
% SDR_n_slow = repmat(n,length(siteList),1);
% SDR_m_slow = squeeze(totalmMatrix(1,1,:))';
% 
% load Seevers_maurer_bestFit_1015_m0_n1.mat
% sites = siteList;
% Seevers_b_slow = squeeze(totalbMatrix(1,1,:))';
% Seevers_n_slow = n;
% Seevers_n_slow = repmat(n,length(siteList),1);
% Seevers_m_slow = squeeze(totalmMatrix(1,1,:))';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% K estimates
SDR_K = @(b,m,n,phi,T2ML) b.*(phi.^m).*(T2ML).^n;
Seevers_K = @(b,m,n,T2ML,T2B,phi) b.*(phi).^m.*((T2ML.^(-1) - T2B.^(-1)).^(-1)).^n;

for kk = 1:length(sites)
    [T2dist, T2logbins, nmrName] = loadRawNMRdata(sites{kk});
   
    [d, K, T2ML, phi, z, SumEch, logK, logT2ML, logPhi, SumEch_3s, SumEch_twm, ...
    SumEch_twm_3s] = loadnmrdata2(nmrName);

    SDR_Kest_slow{kk} = SDR_K(SDR_b_slow(kk),SDR_m_slow(kk),SDR_n_slow(kk),phi,T2ML);
    
    Seevers_Kest_slow{kk} = Seevers_K(Seevers_b_slow(kk),Seevers_m_slow(kk),Seevers_n_slow(kk),...
        T2ML,2.2,phi);
    
    SDR_Kest_fast{kk} = SDR_K(SDR_b_fast(kk),SDR_m_fast(kk),SDR_n_fast(kk),phi,T2ML);
    
    Seevers_Kest_fast{kk} = Seevers_K(Seevers_b_fast(kk),Seevers_m_fast(kk),Seevers_n_fast(kk),...
        T2ML,2.2,phi);
    
    Ktot{kk} = K;
       
end

K_all = vertcat(Ktot{:});
k_estimates = [vertcat(SDR_Kest_fast{:}) vertcat(Seevers_Kest_fast{:}) vertcat(SDR_Kest_slow{:}) vertcat(Seevers_Kest_slow{:})];
k_names = {'SDR','Seevers'};

kDiffFactorTemp = estimateKdiffFactor(K_all, k_estimates, 1);
kDiffFactorMean{:,kk} = mean(kDiffFactorTemp);
kDiffFactor{:,kk} = kDiffFactorTemp;

figure(1)
subplot(2,1,2)
plotKestKdpp_v2(K_all,k_estimates,k_estimates,k_names)
title('Comparing models of K')

title('Kansas + Washington')
%title('Wisconsin')

legend('SDR Fast','Seevers Fast','SDR Slow','Seevers Slow')

ylim([10^-7,10^-2])
xlim([10^-7,10^-2])