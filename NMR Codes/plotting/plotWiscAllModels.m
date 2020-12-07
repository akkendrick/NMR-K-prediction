clear
%close all

load 'wisc_bboot_n2_m0.mat'
load 'Seevers_bestFit_table.mat'
load 'TC_bestFit_240_table.mat'

%sites = {'Site1-WellG5sep','Site1-WellG6sep','Site2-WellPN1','Site2-WellPN2'};
%sites = {'Site1-WellG5','Site1-WellG6'};
sites = {'wisc_all'};   


waterTable = [2.0469,2.1248,5.0285,4.7476]; % rel ground surface
%waterTable = [3.004,2.963,5.727,5.408]; % rel top of casing
depthOffsets = [0.75,0.95,0.75,0.75];

% Sonic depths are relative to ground surface
load('sonicCoreT2B.mat','sonicCoreT2BData')
T2B_depth = sonicCoreT2BData.Depthm;
T2B_depth = flipud(T2B_depth);
 
T2B_peak = sonicCoreT2BData.T2Bpeak;
T2B_peak = flipud(T2B_peak); 

bIndustry = [0.034 0.0425];
mIndustry = [4 4];
nIndustry = [2 2];

bMaurer = [0.80 4.70 0.05 0.11 0.015 0.036 7.3e-4 0.0049];
mMaurer = [4 4 1 1 0 0 0 0];
nMaurer = [2 2 2 2 2 2 1 1];

% bMaurer = [0.015 0.036 7.3e-4 0.0049];
% mMaurer = [0 0 0 0];
% nMaurer = [2 2 1 1];

bKnight = [0.80 5.70];
mKnight = [4 4];
nKnight = [2 2];

bDlubac = [2.4*10^-2];
mDlubac = [2];
nDlubac = [2];

bParsekian = [1.60 2.18 1.18e-1 2.26e-1 5.20e-2 6.84e-2];
mParsekian = [4 4 2 2 1 1];
nParsekian = [2 2 2 2 2 2];

bRen = [1.25e-2 4.04e-2 1.37 2.39e-1 1587.8 86.59...
    1.70e-2 7.00e-4 5.63e-2 1.04e-2 3.78 65.08];
mRen = [1 1 2 2 4 4 1 1 2 2 4 4];
nRen = [2 2 2 2 2 2 1 1 1 1 1 1];

model = 'bootstrap';
%%
    plotIndex = 1;
    figIndex = 1;
SDR_K = @(b,m,n,phi,T2ML) b.*(phi.^m).*(T2ML).^n;
for kk = 1:length(sites)
    site = sites{kk}
    baseDir = '/Volumes/GoogleDrive/My Drive/USGS Project/USGS Data/';
    baseNMRDir = '/Volumes/GoogleDrive/My Drive/USGS Project/NMR-K-prediction/';
 
%     baseDir = 'I:\My Drive\USGS Project\USGS Data\';
%     baseNMRDir = 'I:\My Drive\USGS Project\NMR-K-prediction\';

    T2dist = [];
    T2logbins = [];
    dparam = [];
    kk
    clear Kindustry kTC lK KMat T2MLMat phiMat zMat SumEchMat logKMat logT2MLMat logPhiMat NMRphi kSDR interpT2B kSOE kKGM
    
    if strcmp(site,'wisc_all')
        name = 'wisc_all';
        nmrName = name;
                
        [d, K, T2ML, phi, z, SumEch, logK, logT2ML, logPhi, SumEch_3s, SumEch_twm, ...
        SumEch_twm_3s] = loadnmrdata3(nmrName,baseNMRDir);  
    
        KMat{:,1} = K;
        T2MLMat{:,1} = T2ML;
        phiMat{:,1} = phi;
        zMat{:,1} = z;
        SumEchMat{:,1} = SumEch;
        logKMat{:,1} = logK;
        logT2MLMat{:,1} = logT2ML;
        logPhiMat{:,1} = logPhi;
        

        siteListInd = [1];
        if strcmp(model,'bootstrap')
            SDRb = SDR_Bootstrap(1,siteListInd);
            SDRn = SDR_Bootstrap(2,siteListInd);
            SDRm = SDR_Bootstrap(3,siteListInd);
            
            Seeversb = Seevers_Bootstrap(1, siteListInd);
            Seeversn = Seevers_Bootstrap(2, siteListInd);
            Seeversm = Seevers_Bootstrap(3, siteListInd);
            
            SOEb = [0.0159224451197923];
            SOEn = [1];
            
            KGMrho = [0.0077];
            KGMtau = [2.7227];
            
            TC33b = TC_Bootstrap(1, siteListInd);
            TC33n = TC_Bootstrap(2, siteListInd);
            TC33m = TC_Bootstrap(3, siteListInd);
        elseif strcmp(model,'mcmc')
            SDRb = SDR_MCMC(1,siteListInd);
            SDRn = SDR_MCMC(2,siteListInd);
            SDRm = SDR_MCMC(3,siteListInd);
        else
            SDRb = SDR_Direct(1,siteListInd);
            SDRn = SDR_Direct(2,siteListInd);
            SDRm = SDR_Direct(3,siteListInd);
        end
    end
    
       for j = 1:length(SDRn)
        kSDR{j} = SDR_K(SDRb(j),SDRm(j),SDRn(j),phiMat{j},T2MLMat{j});
     end
     
     finalParameters{1,kk} = SDRb;
     finalParameters{2,kk} = SDRm;
     finalParameters{3,kk} = SDRn;
   
    k_est_n1_ind = [];
    k_est_n2_ind = [];
        
    [k_est_n1_maurer, k_est_n2_maurer, avgK_n1_maurer, avgK_n2_maurer] = computeAvgKModel(siteListInd,bMaurer,mMaurer,nMaurer,phiMat,T2MLMat);
    [k_est_n1_dlubac, k_est_n2_dlubac, avgK_n1_dlubac, avgK_n2_dlubac] = computeAvgKModel(siteListInd,bDlubac,mDlubac,nDlubac,phiMat,T2MLMat);

    k_estimates_best = [k_est_n1_maurer k_est_n2_dlubac];
   % k_avg_best = [avgK_n1_maurer avgK_n2_dlubac vertcat(kSDR{:})];
    
    k_avg_best = [vertcat(kSDR{:})];
    k_names = [{'Best SDR'}]
    
        colors = {[0.9290,0.6940,0.1250],[0.466,0.6740,0.1880]};
    figure(2)
    plotKestKdpp_wdepth(vertcat(KMat{:}),k_estimates_best,k_avg_best,vertcat(zMat{:}),k_names,colors)   
    title('Estimating K with all Wisconsin data')  
    figIndex = figIndex + 1;
    set(gca,'FontSize',14)
    
end