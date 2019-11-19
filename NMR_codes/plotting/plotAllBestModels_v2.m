% Plot all best models, improved 
% Plot results of prior estimates of bounding values of b,m,n

clear
%close all

SDR = load('SDR_bestFit_1014_m0_n2.mat');
%Seevers = load('Seevers_bestFit_1101_m0_n2.mat');
Seevers = load('Seevers_model_pairs_noSep_m0_n2.mat')

load 'TC_bestFit_240_table.mat'

sites = {'Site1-WellG5','Site1-WellG6','Site2-WellPN1','Site2-WellPN2'};

waterTable = [2.0469,2.1248,5.0285,4.7476]; % rel ground surface
%waterTable = [3.004,2.963,5.727,5.408]; % rel top of casing
depthOffsets = [0.75,0.95,0.75,0.75];

% Sonic depths are relative to ground surface
load('sonicCoreT2B.mat','sonicCoreT2BData')
T2B_depth = sonicCoreT2BData.Depthm;
T2B_depth = flipud(T2B_depth);
 
T2B_peak = sonicCoreT2BData.T2Bpeak;
T2B_peak = flipud(T2B_peak); 
goodT2B = T2B_peak(T2B_peak > 1000);
T2Bavg = mean(goodT2B)*10^-3;

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
 TC_K = @(c,m,n,phi,indexQuotient) c.*(phi).^m.*(indexQuotient).^n;
 SOE_K = @(b,n,SOE) b.*(SOE).^n;
 Seevers_K = @(b,m,n,T2ML,T2B,phi) b.*(phi).^m.*((T2ML.^(-1) - T2B.^(-1)).^(-1)).^n;
 
  
Temp = 20;  % temperature in degress C 
rho = @(Tt) 1000*(1 - ((Tt+288.94)./(508929*(Tt+68.12))).*(Tt-3.98).^2); % kg/m^3
eta = @(Tt) 0.0013 - 1.7e-5*Tt;         % Pa -s
%Tb = @(Tt) 3.3 + 0.044*(Tt - 35);       % seconds

Tb = @(Tt) 200;
D = @(Tt) (1.0413 + 0.039828*Tt + 0.00040318*Tt.^2).*1e-9;  % m^2/s 
g = 9.8;    %m/s^2
tort = 1/(1.5^2); 
t1 = (rho(Temp)*g)/(8*eta(Temp)); % 

num2 = @(T2) 4*D(Temp)*Tb(Temp)*T2;
denom2 = @(T2) Tb(Temp) - T2; 
    
f12 = @(rho) (D(Temp)./rho);  
SQterm = @(rho,T2) sqrt(f12(rho).^2 + (num2(T2)./denom2(T2))); 

KGM_lK = @(rho,tau,m,lphi,T2) log10(1/tau^2) + log10(t1) + m*lphi + 2*log10(SQterm(rho,T2)-f12(rho)); 

 
for kk = 1:length(sites)
    site = sites{kk};
     baseDir = '/Volumes/GoogleDrive/My Drive/Stanford/USGS Project/Field Data/USGS Data/';
     baseNMRDir = '/Volumes/GoogleDrive/My Drive/Stanford/USGS Project/Field Data/NMR-K-prediction/';
%  
%    baseDir = 'I:\My Drive\USGS Project\USGS Data\';
%    baseNMRDir = 'I:\My Drive\USGS Project\NMR-K-prediction\';

    T2dist = [];
    T2logbins = [];
    dparam = [];
    kk
    clear Kindustry kTC lK KMat T2MLMat phiMat zMat SumEchMat logKMat logT2MLMat logPhiMat NMRphi kSDR interpT2B kSOE kKGM
    
    n_global = 2;
    m_global = 0;
    
    if strcmp(site,'Site1-WellG5')
        name = 'G5_W1_tr5_20x_16p5_up_F1n2_wRIN_wRFI_Reg50_Va1';
        nmrName = name;

        in1 = [baseDir site '/' name '/' name '_T2_dist' '.txt']; 
        in2 = [baseDir site '/' name '/' name '_T2_bins_log10s' '.txt']; 
        in3 = [baseDir site '/' name '/' name '_1Dvectors' '.txt'];

        T2dist = load(in1);
        T2dist(:,1) = T2dist(:,1) - depthOffsets(1);
        T2logbins = load(in2);
        
        dparam = dlmread(in3,'\t',1,0); 


        [d, K, T2ML, phi, z, SumEch, logK, logT2ML, logPhi, SumEch_3s, SumEch_twm, ...
        SumEch_twm_3s] = loadnmrdata2(nmrName); 

        KMat{:,1} = K;
        T2MLMat{:,1} = T2ML;
        phiMat{:,1} = phi;
        zMat{:,1} = z;
        SumEchMat{:,1} = SumEch;
        logKMat{:,1} = logK;
        logT2MLMat{:,1} = logT2ML;
        logPhiMat{:,1} = logPhi;

        NMRphi(:,1) = dparam(:,2);

        indexQuotient{:,1} = TC_indexQuotient{4}; 
        
        nmrDepths = z(1:end,1);
        interpT2B = interp1(T2B_depth, T2B_peak, nmrDepths)*10^-3; %Convert from ms to s
        
        
        
        siteListInd = 2;
        if strcmp(model,'bootstrap')
            SDRb = SDR.totalbMatrix(1,1,siteListInd);
            SDRn = n_global;
            SDRm = m_global; 
            
            Seeversb = Seevers.totalbMatrix(1,1,siteListInd);
            Seeversn = n_global;
            Seeversm = m_global;

            SOEb = [0.00512209505212264];
            SOEn = 1;
            
            KGMtau = 1.7378;
            KGMrho = 5.533501092157371e-05;
            
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
    
        
    elseif strcmp(site,'Site1-WellG6')
        name = 'G6_W2_tr5_20x_16p75_up_F_wRIN_wRFI_reg50_Va1';
        nmrName = name;

        in1 = [baseDir site '/' name '/' name '_T2_dist' '.txt']; 
        in2 = [baseDir site '/' name '/' name '_T2_bins_log10s' '.txt']; 
        in3 = [baseDir site '/' name '/' name '_1Dvectors' '.txt'];

        T2dist = load(in1);
        T2dist(:,1) = T2dist(:,1) - depthOffsets(2);
        
        dparam = dlmread(in3,'\t',1,0); 
        T2logbins = load(in2);


        [d, K, T2ML, phi, z, SumEch, logK, logT2ML, logPhi, SumEch_3s, SumEch_twm, ...
        SumEch_twm_3s] = loadnmrdata2(nmrName);  

        KMat{:,1} = K;
        T2MLMat{:,1} = T2ML;
        phiMat{:,1} = phi;
        zMat{:,1} = z;
        SumEchMat{:,1} = SumEch;
        logKMat{:,1} = logK;
        logT2MLMat{:,1} = logT2ML;
        logPhiMat{:,1} = logPhi;

        NMRphi(:,1) = dparam(:,2);

        indexQuotient{:,1} = TC_indexQuotient{1}; 

        
        siteListInd = [1];
        if strcmp(model,'bootstrap')
            SDRb = SDR.totalbMatrix(1,1,siteListInd);
            SDRn = n_global;
            SDRm = m_global; 
            
            Seeversb = Seevers.totalbMatrix(1,1,siteListInd);
            Seeversn = n_global;
            Seeversm = m_global;

            
            SOEb = [0.00451558405202521];
            SOEn = [1];
            
            KGMrho = 2.060629913270001e-05;
            KGMtau = 1;
            
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
        nmrDepths = vertcat(zMat{:});
        interpT2B = interp1(T2B_depth, T2B_peak, nmrDepths)*10^-3; %Convert from ms to s
        
    elseif strcmp(site,'Site2-WellPN1')
        name = 'Pl_W1_Tr5_20x_MPp75aLS_F1n2_wRIN_wRFI_Reg50_Va1';
        nmrName = name;

        interpT2B = [];
        
        in1 = [baseDir site '/' name '/' name '_T2_dist' '.txt']; 
        in2 = [baseDir site '/' name '/' name '_T2_bins_log10s' '.txt']; 
        in3 = [baseDir site '/' name '/' name '_1Dvectors' '.txt'];

        T2dist = load(in1);
        T2dist(:,1) = T2dist(:,1) - depthOffsets(3);
        T2logbins = load(in2);

        dparam = dlmread(in3,'\t',1,0); 


        [d, K, T2ML, phi, z, SumEch, logK, logT2ML, logPhi, SumEch_3s, SumEch_twm, ...
        SumEch_twm_3s] = loadnmrdata2(nmrName);  

        KMat{:,1} = K;
        T2MLMat{:,1} = T2ML;
        phiMat{:,1} = phi;
        zMat{:,1} = z;
        SumEchMat{:,1} = SumEch;
        logKMat{:,1} = logK;
        logT2MLMat{:,1} = logT2ML;
        logPhiMat{:,1} = logPhi;

        NMRphi(:,1) = dparam(:,2);
        
        indexQuotient{:,1} = TC_indexQuotient{7}; 

        siteListInd = [3];
        %siteListInd = [7];
        if strcmp(model,'bootstrap')
            SDRb = SDR.totalbMatrix(1,1,siteListInd);
            SDRn = n_global;
            SDRm = m_global; 
            
            Seeversb = Seevers.totalbMatrix(1,1,siteListInd);
            Seeversn = n_global;
            Seeversm = m_global;

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

    elseif strcmp(site,'Site2-WellPN2')
        name = 'W2_Tr5_20x_MPp75aLS_Reg50_wRIN_wRFI_Va1';
        nmrName = name;

        in1 = [baseDir site '/' name '/' name '_T2_dist' '.txt']; 
        in2 = [baseDir site '/' name '/' name '_T2_bins_log10s' '.txt']; 
        in3 = [baseDir site '/' name '/' name '_1Dvectors' '.txt'];

        interpT2B = [];
        T2logbins = load(in2);

        T2dist = load(in1);
        T2dist(:,1) = T2dist(:,1) - depthOffsets(4);
        
        dparam = dlmread(in3,'\t',1,0); 

        [d, K, T2ML, phi, z, SumEch, logK, logT2ML, logPhi, SumEch_3s, SumEch_twm, ...
        SumEch_twm_3s] = loadnmrdata2(nmrName);  

        KMat{:,1} = K;
        T2MLMat{:,1} = T2ML;
        phiMat{:,1} = phi;
        zMat{:,1} = z;
        SumEchMat{:,1} = SumEch;
        logKMat{:,1} = logK;
        logT2MLMat{:,1} = logT2ML;
        logPhiMat{:,1} = logPhi;

        NMRphi(:,1) = dparam(:,2);
        indexQuotient{:,1} = TC_indexQuotient{8}; 
        %siteListInd = [8];
        siteListInd = [4];
        
        if strcmp(model,'bootstrap')
            SDRb = SDR.totalbMatrix(1,1,siteListInd);
            SDRn = n_global;
            SDRm = m_global; 
            
            Seeversb = Seevers.totalbMatrix(1,1,siteListInd);
            Seeversn = n_global;
            Seeversm = m_global;

            SOEb = [0.009135553542862];
            SOEn = [1];
            
            KGMtau = [3.467];
            KGMrho = 1;
            
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

%%

    noPorosity = 1;  

%%

    k_est_n1_ind = [];
    k_est_n2_ind = [];
        
    TC_K = @(c,m,n,phi,indexQuotient) c.*(phi).^m.*(indexQuotient).^n;
 SOE_K = @(b,n,SOE) b.*(SOE).^n;
    
     for j = 1:length(SDRn)
        kSDR{j} = SDR_K(SDRb(j),SDRm(j),SDRn(j),phiMat{j},T2MLMat{j});
        kSeevers{j} = Seevers_K(Seeversb(j),Seeversm(j),Seeversn(j),T2MLMat{j},T2Bavg,phiMat{j});
        kSOE{j} = SOE_K(SOEb(j),SOEn(j),SumEchMat{j});
        kTC{j} = TC_K(TC33b(j),TC33m(j),TC33n(j),phiMat{j},indexQuotient{:,j});
        
        lkKGM{j} = KGM_lK(KGMrho(j), KGMtau(j), 0, logPhiMat{j}, T2MLMat{j}); 
        kKGM{j} = 10.^lkKGM{j};
     end
     
     finalParameters{1,kk} = SDRb;
     finalParameters{2,kk} = SDRm;
     finalParameters{3,kk} = SDRn;

    zTest = vertcat(zMat{:});
    [z, index] = sort(vertcat(zMat{:}),'descend');

    % SORT EVERYTHING
    K = vertcat(KMat{:});
    K = K(index);
    
    %error_all(kk,:) = computeError(K,k_avg_all);

    kSDR_temp = vertcat(kSDR{:});
    kSeevers_temp = vertcat(kSeevers{:});
    kKGM_temp = vertcat(kKGM{:});
    kTC_temp = vertcat(kTC{:});
    kSOE_temp = vertcat(kSOE{:});
    
    bestSDRModel{kk} = kSDR_temp(index);
    bestSeeversModel{kk} = kSeevers_temp(index);
    bestKGMModel{kk} = kKGM_temp(index);
    bestTCModel{kk} = kTC_temp(index);
    bestSOEModel{kk} = kSOE_temp(index);
    
    SDRError{:,kk} = computeError(K,bestSDRModel{kk});
    SDRErrorFactor{:,kk} = mean(estimateKdiffFactor(K,bestSDRModel{kk},1));
    
    SeeversError{:,kk} = computeError(K,bestSeeversModel{kk});
    SeeversErrorFactor{:,kk} = mean(estimateKdiffFactor(K,bestSeeversModel{kk},1));
    
    KGMError{:,kk} = computeError(K,bestKGMModel{kk});
    KGMErrorFactor{:,kk} = mean(estimateKdiffFactor(K,bestKGMModel{kk},1));
    
    TCError{:,kk} = computeError(K,bestTCModel{kk});
    TCErrorFactor{:,kk} = mean(estimateKdiffFactor(K,bestTCModel{kk},1));
    
    SOEError{:,kk} = computeError(K,bestSOEModel{kk});
    SOEErrorFactor{:,kk} = mean(estimateKdiffFactor(K,bestSOEModel{kk},1));
    
    k_est_all = [bestSDRModel{kk},bestSeeversModel{kk},bestKGMModel{kk},...
        bestTCModel{kk},bestSOEModel{kk}];
    k_avg_all = [bestSDRModel{kk},bestSeeversModel{kk},bestKGMModel{kk},...
        bestTCModel{kk},bestSOEModel{kk}];
    k_names_all = [{'SDR'},{'Seevers'},{'KGM'},{'TC'},{'SOE'}];
    
%    % plotKestKdpp_v2(vertcat(KMat{:}),k_estimates_n2,k_avg_n2,k_names_n2)   
%     %title(site)

    figure(20)   
    hold on
    subplot(4,1,kk)
    plotKestKdpp_wdepth(vertcat(KMat{:}),k_est_all,k_avg_all,vertcat(zMat{:}),k_names_all)   
    title(site) 
    
    figure(21)
    
    plotKPhiwithDepth_v2(vertcat(KMat{:}),vertcat(phiMat{:}),waterTable(kk),vertcat(zMat{:}),...
        T2dist,T2logbins,k_est_all,k_avg_all,k_names_all,1,kk)
%     

    legend(k_names_all)

%     figIndex = figIndex + 1;
%     
%     figure(2)
%     subplot(4,3,figIndex)
%     plotKFactor_wdepth(vertcat(KMat{:}),k_estimates_n2,k_avg_n2,vertcat(zMat{:}),k_names_n2)   
%     title(site) 
% 
%     figIndex = figIndex + 2;
    

%     subplot(2,5,5.5);
%     %poshL = get(hL,'position'); 
%     legend(k_names_n2);
%     set(gca,'FontSize',16)
%     %set(lgd,'position',poshL);      % Adjusting legend's position
%     %axis(hL,'off');                 % Turning its axis off
%     axis off
%     
%    
%     plotKPhiwithDepth_v2(vertcat(KMat{:}),NMRphi,waterTable(kk),vertcat(zMat{:}),T2dist,...
%     T2logbins,k_estimates_n1,k_avg_n1,k_names_n1,noPorosity)
%     title(site)
%     
%     plotKestKdpp_v2(vertcat(KMat{:}),k_estimates_n1,k_avg_n1,k_names_n1)   
%     title(site) 
%     
%     plotKFactor_v2(vertcat(KMat{:}),k_estimates_n1,k_avg_n1,k_names_n1)   
    title(site)  
    



end

figure(22)

box on
grid on
hold on
scatter([1:4],vertcat(SDRErrorFactor{:}),40,'filled')
scatter([1:4],vertcat(SeeversErrorFactor{:}),40,'filled')
scatter([1:4],vertcat(KGMErrorFactor{:}),40,'filled')
scatter([1:4],vertcat(TCErrorFactor{:}),40,'filled')
scatter([1:4],vertcat(SOEErrorFactor{:}),40,'filled')

ylim([0 10])
set(gca,'yscale','log')
legend(k_names_all)
%save('bestKModels.mat','dlubacModel','bestSDRModel')

%legend(k_names_n2);
