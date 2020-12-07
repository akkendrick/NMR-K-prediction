% Plot results of prior estimates of bounding values of b,m,n

clear
%close all

%load 'SDR_bestFit_table_m0_n2.mat'
%load 'Seevers_bestFit_table.mat'
%load 'TC_bestFit_240_table.mat'

load SDR_bestFit_1122_m0_n2.mat
%load Seevers_bestFit_1201_m1_n2_T2BAvg.mat

%sites = {'Site1-WellG5sep','Site1-WellG6sep','Site2-WellPN1','Site2-WellPN2'};
%sites = {'Site1-WellG5','Site1-WellG6'};
sites = {'Site1-WellG5','Site1-WellG6','Site2-WellPN1','Site2-WellPN2'};
names = {'Adams Well G5','Adams Well G6','Plainfield Lake Well PN1','Plainfield Lake Well PN2'};
%sites = {'wisc_all'};   

waterTable = [1.181,1.181,5.0285,4.7476];
%waterTable = [2.0469,2.1248,5.0285,4.7476]; % rel ground surface
%waterTable = [3.004,2.963,5.727,5.408]; % rel top of casing
depthOffsets = [0.75,0.95,0.75,0.75];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prior model parameters 
SDR_b = [0.0040 0.0024 0.0075 0.0047];
SDR_m = [0 0 0 0];
SDR_n = [2 2 2 2];

% Seevers
Seevers_b = [0.0016 0.0032 0.0058 0.0038];
Seevers_n = [2 2 2 2];
Seevers_m = [0 0 0 0];

% KGM (has been checked, looks good)
KGM_tau = [1 1.7378 2.7227 3.4674];
KGM_rho = [2.0606e-05 5.5335e-05 0.0077 1];
KGM_m = [0 0 0 0];

% SOE
SOE_n = [1 1 1 1];
SOE_b = [0.0045 0.0052 0.0158 0.0092];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% K Models
SDR_K = @(b,m,n,phi,T2ML) b.*(phi.^m).*(T2ML).^n;
SOE_K = @(b,n,SOE) b.*(SOE).^n;
Seevers_K = @(b,m,n,T2ML,T2B,phi) b.*(phi).^m.*((T2ML.^(-1) - T2B.^(-1)).^(-1)).^n;

Temp = 20;  % temperature in degress C 
rho = @(Tt) 1000*(1 - ((Tt+288.94)./(508929*(Tt+68.12))).*(Tt-3.98).^2); % kg/m^3
eta = @(Tt) 0.0013 - 1.7e-5*Tt;         % Pa -s
Tb = @(Tt) 3.3 + 0.044*(Tt - 35);       % seconds
D = @(Tt) (1.0413 + 0.039828*Tt + 0.00040318*Tt.^2).*1e-9;  % m^2/s 
g = 9.8;    %m/s^2
tort = 1/(1.5^2); 
t1 = (rho(Temp)*g)/(8*eta(Temp)); % 

num2 = @(T2) 4*D(Temp)*Tb(Temp)*T2;
denom2 = @(T2) Tb(Temp) - T2; 
    
f12 = @(rho) (D(Temp)./rho);  
SQterm = @(rho,T2) sqrt(f12(rho).^2 + (num2(T2)./denom2(T2))); 

KGM_lK = @(rho,tau,m,lphi,T2) log10(1/tau^2) + log10(t1) + m*lphi + 2*log10(SQterm(rho,T2)-f12(rho)); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sonic depths are relative to ground surface
load('sonicCoreT2B.mat','sonicCoreT2BData')
T2B_depth = sonicCoreT2BData.Depthm;
T2B_depth = flipud(T2B_depth);
 
T2B_peak = sonicCoreT2BData.T2Bpeak;
T2B_peak = flipud(T2B_peak); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bIndustry = [0.034];
mIndustry = [4];
nIndustry = [2];

bMaurer = [0.022 0.036 0.027 0.016 0.024 0.032 0.028];
mMaurer = [0 0 0 0 0 0 0];
nMaurer = [2 2 2 2 2 2 2];

bDlubac = [2.4*10^-2];
mDlubac = [2];
nDlubac = [2];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
model = 'bootstrap';
%%
    plotIndex = 1;
    figIndex = 1;
SDR_K = @(b,m,n,phi,T2ML) b.*(phi.^m).*(T2ML).^n;
for kk = 1:length(sites)
    site = sites{kk}
    baseDir = '/Volumes/GoogleDrive/My Drive/Stanford/USGS Project/Field Data/USGS Data/';
    baseNMRDir = '/Volumes/GoogleDrive/My Drive/Stanford/USGS Project/NMR-K-prediction/';
%  
%     baseDir = 'I:\My Drive\USGS Project\USGS Data\';
%     baseNMRDir = 'I:\My Drive\USGS Project\NMR-K-prediction\';

    T2dist = [];
    T2logbins = [];
    dparam = [];
    kk
    clear Kindustry kTC lK KMat T2MLMat phiMat zMat SumEchMat logKMat logT2MLMat logPhiMat NMRphi kSDR interpT2B kSOE kKGM
    
    if strcmp(site,'Site1-WellG5')
        name = 'G5_W1_tr5_20x_16p5_up_F1n2_wRIN_wRFI_Reg50_Va1';
        nmrName = name;

        in1 = [baseDir site '/' name '/' name '_T2_dist' '.txt']; 
        in2 = [baseDir site '/' name '/' name '_T2_bins_log10s' '.txt']; 
        in3 = [baseDir site '/' name '/' name '_1Dvectors' '.txt'];

        T2dist = load(in1) - depthOffsets(1); 
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

       % indexQuotient{:,1} = TC_indexQuotient{4}; 
        
        nmrDepths = z(1:end,1);
        interpT2B = interp1(T2B_depth, T2B_peak, nmrDepths)*10^-3; %Convert from ms to s
        
        siteListInd = 2;
        %siteListInd = 4;
        if strcmp(model,'bootstrap')
            SDRb = SDR_Bootstrap(1,siteListInd);
            SDRn = SDR_Bootstrap(2,siteListInd);
            SDRm = SDR_Bootstrap(3,siteListInd);
            
            Seeversb = Seevers_Bootstrap(1, siteListInd);
            Seeversn = Seevers_Bootstrap(2, siteListInd);
            Seeversm = Seevers_Bootstrap(3, siteListInd);

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
    elseif  strcmp(site,'Site1-WellG5sep')
        site = 'Site1-WellG5';
        name1 = 'G5_W1_tr5_20x_16p5_up_F1n2_wRIN_wRFI_Reg50_Va1';
        
        nmrName1 = 'G5_W1_tr5_20x_16p5_up_F1n2_wRIN_wRFI_Reg50_Va1_above';

        [d, K, T2ML, phi, z, SumEch, logK, logT2ML, logPhi, SumEch_3s, SumEch_twm, ...
        SumEch_twm_3s] = loadnmrdata3(nmrName1,baseNMRDir); 

        KMat{:,1} = K;
        T2MLMat{:,1} = T2ML;
        phiMat{:,1} = phi;
        zMat{:,1} = z;
        SumEchMat{:,1} = SumEch;
        logKMat{:,1} = logK;
        logT2MLMat{:,1} = logT2ML;
        logPhiMat{:,1} = logPhi;

        name2 = 'G5_W1_tr5_20x_16p5_up_F1n2_wRIN_wRFI_Reg50_Va1';

        in1 = [baseDir site '/' name2 '/' name2 '_T2_dist' '.txt']; 
        in2 = [baseDir site '/' name2 '/' name2 '_T2_bins_log10s' '.txt']; 
        in3 = [baseDir site '/' name2 '/' name2 '_1Dvectors' '.txt'];

        T2dist = load(in1) - depthOffsets(1);
        T2logbins = load(in2);
        dparam = dlmread(in3,'\t',1,0); 

        nmrName2 = 'G5_W1_tr5_20x_16p5_up_F1n2_wRIN_wRFI_Reg50_Va1_below';


        [d, K, T2ML, phi, z, SumEch, logK, logT2ML, logPhi, SumEch_3s, SumEch_twm, ...
        SumEch_twm_3s] = loadnmrdata3(nmrName2,baseNMRDir); 

        KMat{:,2} = K;
        T2MLMat{:,2} = T2ML;
        phiMat{:,2} = phi;
        zMat{:,2} = z;
        SumEchMat{:,2} = SumEch;
        logKMat{:,2} = logK;
        logT2MLMat{:,2} = logT2ML;
        logPhiMat{:,2} = logPhi;

        NMRphi = dparam(:,2);
        
        %indexQuotient{:,1} = TC_indexQuotient{5}; 
        %indexQuotient{:,2} = TC_indexQuotient{6}; 

        
        %siteListInd = [5 6];
        if strcmp(model,'bootstrap')
            SDRb = SDR_Bootstrap(1,siteListInd);
            SDRn = SDR_Bootstrap(2,siteListInd);
            SDRm = SDR_Bootstrap(3,siteListInd);

            Seeversb = Seevers_Bootstrap(1, siteListInd);
            Seeversn = Seevers_Bootstrap(2, siteListInd);
            Seeversm = Seevers_Bootstrap(3, siteListInd);

            SOEb = [0.00426495739362032 0.0115188123981797];
            SOEn = [1 1];
            
            KGMrho = [8.729713683881130e-05 1];
            KGMtau = [2.483 2.917];
            
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
        
        interpT2B{:,1} = interp1(T2B_depth, T2B_peak, zMat{:,1})*10^-3; %Convert from ms to s
        interpT2B{:,2} = interp1(T2B_depth, T2B_peak, zMat{:,2})*10^-3; %Convert from ms to s
        
    elseif strcmp(site,'Site1-WellG6')
        name = 'G6_W2_tr5_20x_16p75_up_F_wRIN_wRFI_reg50_Va1';
        nmrName = name;

        in1 = [baseDir site '/' name '/' name '_T2_dist' '.txt']; 
        in2 = [baseDir site '/' name '/' name '_T2_bins_log10s' '.txt']; 
        in3 = [baseDir site '/' name '/' name '_1Dvectors' '.txt'];

        T2dist = load(in1) - depthOffsets(2); 
        T2logbins = load(in2);
        dparam = dlmread(in3,'\t',1,0); 


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

        NMRphi(:,1) = dparam(:,2);

        %indexQuotient{:,1} = TC_indexQuotient{1}; 

        
%        siteListInd = [1];
        siteListInd = 1;

        if strcmp(model,'bootstrap')
            SDRb = SDR_Bootstrap(1,siteListInd)
            SDRn = SDR_Bootstrap(2,siteListInd);
            SDRm = SDR_Bootstrap(3,siteListInd);
           
            Seeversb = Seevers_Bootstrap(1, siteListInd);
            Seeversn = Seevers_Bootstrap(2, siteListInd);
            Seeversm = Seevers_Bootstrap(3, siteListInd);
            
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
        
    elseif strcmp(site,'Site1-WellG6sep')
        site = 'Site1-WellG6';
        name1 = 'G6_W2_tr5_20x_16p75_up_F_wRIN_wRFI_reg50_Va1';
        nmrName1 = 'G6_W2_tr5_20x_16p75_up_F_wRIN_wRFI_reg50_Va1_above';

        [d, K, T2ML, phi, z, SumEch, logK, logT2ML, logPhi, SumEch_3s, SumEch_twm, ...
        SumEch_twm_3s] = loadnmrdata3(nmrName1,baseNMRDir); 

        KMat{:,1} = K;
        T2MLMat{:,1} = T2ML;
        phiMat{:,1} = phi;
        zMat{:,1} = z;
        SumEchMat{:,1} = SumEch;
        logKMat{:,1} = logK;
        logT2MLMat{:,1} = logT2ML;
        logPhiMat{:,1} = logPhi;
       
        name2 = 'G6_W2_tr5_20x_16p75_up_F_wRIN_wRFI_reg50_Va1';

        in1 = [baseDir site '/' name2 '/' name2 '_T2_dist' '.txt']; 
        in2 = [baseDir site '/' name2 '/' name2 '_T2_bins_log10s' '.txt']; 
        in3 = [baseDir site '/' name2 '/' name2 '_1Dvectors' '.txt'];

        T2dist = load(in1) - depthOffsets(2); 
        T2logbins = load(in2);
        dparam = dlmread(in3,'\t',1,0); 

        nmrName2 = 'G6_W2_tr5_20x_16p75_up_F_wRIN_wRFI_reg50_Va1_below';

        [d, K, T2ML, phi, z, SumEch, logK, logT2ML, logPhi, SumEch_3s, SumEch_twm, ...
        SumEch_twm_3s] = loadnmrdata3(nmrName2,baseNMRDir);  

        KMat{:,2} = K;
        T2MLMat{:,2} = T2ML;
        phiMat{:,2} = phi;
        zMat{:,2} = z;
        SumEchMat{:,2} = SumEch;
        logKMat{:,2} = logK;
        logT2MLMat{:,2} = logT2ML;
        logPhiMat{:,2} = logPhi;

        NMRphi = dparam(:,2);

        %indexQuotient{:,1} = TC_indexQuotient{2}; 
        %indexQuotient{:,2} = TC_indexQuotient{3}; 

    %    siteListInd = [2 3];
        if strcmp(model,'bootstrap')
            SDRb = SDR_Bootstrap(1,siteListInd);
            SDRn = SDR_Bootstrap(2,siteListInd);
            SDRm = SDR_Bootstrap(3,siteListInd);
        
            Seeversb = Seevers_Bootstrap(1, siteListInd);
            Seeversn = Seevers_Bootstrap(2, siteListInd);
            Seeversm = Seevers_Bootstrap(3, siteListInd);
            
            SOEb = [0.00336947073635038 0.012235616139267]; 
            SOEn = [1 1];
            
            KGMrho = [1.770108958317422e-05 1];
            KGMtau = [1 2.754];
            
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

        interpT2B{:,1} = interp1(T2B_depth, T2B_peak, zMat{:,1})*10^-3; %Convert from ms to s
        interpT2B{:,2} = interp1(T2B_depth, T2B_peak, zMat{:,2})*10^-3; %Convert from ms to s
        
    elseif strcmp(site,'Site2-WellPN1')
        name = 'Pl_W1_Tr5_20x_MPp75aLS_F1n2_wRIN_wRFI_Reg50_Va1';
        nmrName = name;

        interpT2B = [];
        
        in1 = [baseDir site '/' name '/' name '_T2_dist' '.txt']; 
        in2 = [baseDir site '/' name '/' name '_T2_bins_log10s' '.txt']; 
        in3 = [baseDir site '/' name '/' name '_1Dvectors' '.txt'];

        T2dist = load(in1) - depthOffsets(3); 
        T2logbins = load(in2);
        dparam = dlmread(in3,'\t',1,0); 


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

        NMRphi(:,1) = dparam(:,2);
        
        %indexQuotient{:,1} = TC_indexQuotient{7}; 


      %  siteListInd = [7];
        siteListInd = [3];
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

    elseif strcmp(site,'Site2-WellPN2')
        name = 'W2_Tr5_20x_MPp75aLS_Reg50_wRIN_wRFI_Va1';
        nmrName = name;

        in1 = [baseDir site '/' name '/' name '_T2_dist' '.txt']; 
        in2 = [baseDir site '/' name '/' name '_T2_bins_log10s' '.txt']; 
        in3 = [baseDir site '/' name '/' name '_1Dvectors' '.txt'];

        interpT2B = [];
        
        T2dist = load(in1) - depthOffsets(4); 
        T2logbins = load(in2);
        dparam = dlmread(in3,'\t',1,0); 

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

        NMRphi(:,1) = dparam(:,2);
        %indexQuotient{:,1} = TC_indexQuotient{8}; 
        %siteListInd = [8];
        siteListInd = [4];

        if strcmp(model,'bootstrap')
            SDRb = SDR_Bootstrap(1,siteListInd);
            SDRn = SDR_Bootstrap(2,siteListInd);
            SDRm = SDR_Bootstrap(3,siteListInd);
            
            Seeversb = Seevers_Bootstrap(1, siteListInd);
            Seeversn = Seevers_Bootstrap(2, siteListInd);
            Seeversm = Seevers_Bootstrap(3, siteListInd);

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
    elseif strcmp(site,'wisc_all')
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

%%

    noPorosity = 1;


%% Only plot best model 
    
    
%      for j = 1:length(SDRn)
%         kSDR{j} = SDR_K(SDRb(j),SDRm(j),SDRn(j),phiMat{j},T2MLMat{j});
%      end
     
    kSDR{kk} = SDR_K(SDR_b(kk), SDR_m(kk), SDR_n(kk), phi,T2ML);

     kSeevers{kk} = Seevers_K(Seevers_b(kk),Seevers_m(kk),Seevers_n(kk),...
        T2ML,2.2,phi); 
    
     
     kSOE{kk} = SOE_K(SOE_b(kk),SOE_n(kk),SumEch);
     
     lkKGM{kk} = KGM_lK(KGM_rho(kk),KGM_tau(kk),KGM_m(kk),log10(phi),T2ML);
     kKGM{kk} = 10.^lkKGM{kk};
     
     
    [k_est_n1_ind, k_est_n2_ind, avgK_n1_ind, avgK_n2_ind] = computeAvgKModel(siteListInd,bIndustry,mIndustry,nIndustry,phiMat,T2MLMat);
    [k_est_n1_maurer, k_est_n2_maurer, avgK_n1_maurer, avgK_n2_maurer] = computeAvgKModel(siteListInd,bMaurer,mMaurer,nMaurer,phiMat,T2MLMat);
    [k_est_n1_dlubac, k_est_n2_dlubac, avgK_n1_dlubac, avgK_n2_dlubac] = computeAvgKModel(siteListInd,bDlubac,mDlubac,nDlubac,phiMat,T2MLMat);

    
    
    
    k_estimates_best = [k_est_n2_maurer k_est_n2_dlubac];
     
    finalParameters{1,kk} = SDRb;
    finalParameters{2,kk} = SDRm;
    finalParameters{3,kk} = SDRn;
   
    k_est_n1_ind = [];
    k_est_n2_ind = [];
        
   % k_avg_best = [avgK_n1_maurer avgK_n2_dlubac vertcat(kSDR{:})];
    k_avg_best = [vertcat(kSDR{:})];% kSeevers{kk} kSOE{kk} kKGM{kk}];
    k_names = [{'Best SDR'}]
    
   % k_names = [{'Maurer'} {'Dlubac'} {'Best SDR'}];
    %k_sym = [{'.'} {'.'} {'*'} {'^'} {'^'} {'^'} {'^'} {'^'} {'^'} {'^'} {'>'}];

    plotKPhiwithDepth_v2(vertcat(KMat{:}),NMRphi,waterTable(kk),vertcat(zMat{:}),T2dist,...
    T2logbins,k_estimates_best,k_avg_best,k_names,noPorosity,plotIndex)
    title(names{kk})
    
    fileName = strcat(sites{kk},'_SDRprofile.svg');
    print('-dsvg','-r300',fileName)
    
    plotIndex = plotIndex + 2;

%     plotKestKdpp_depth(vertcat(KMat{:}),k_estimates_n2,k_avg_n2,k_names)   
%     title(site) 
%     
%     plotKFactor_depth(vertcat(KMat{:}),k_estimates_n2,k_avg_n2,k_names)   
%     title(site)     

    %k_avg_best = [vertcat(kSDR{:}) kSeevers{kk} kSOE{kk} kKGM{kk}];
    
    %k_avg_best = [vertcat(kSDR{:})  avgK_n2_dlubac avgK_n2_ind avgK_n2_maurer];
    k_avg_best = [vertcat(kSDR{:}) avgK_n2_maurer avgK_n2_dlubac avgK_n2_ind];

    figure(30)
    subplot(4,1,figIndex)
    plotKestKdpp_wdepth(vertcat(KMat{:}),k_estimates_best,k_avg_best,vertcat(zMat{:}),k_names)   
    title(site)  
    figIndex = figIndex + 1;
    
%     subplot(4,2,figIndex)
%     plotKFactor_wdepth(vertcat(KMat{:}),k_estimates_best,k_avg_best,vertcat(zMat{:}),k_names)   
%     figIndex = figIndex + 1;
%     title(site)
    
    kDiffFactorAvg = mean(estimateKdiffFactor(vertcat(KMat{:}),avgK_n2_maurer,0));
    
    kDiffFactor = (estimateKdiffFactor(vertcat(KMat{:}),k_avg_best,1));
    mean(kDiffFactor)
    
    
    k_Maurer_all{kk} = avgK_n2_maurer;
    k_Dlubac_all{kk} = avgK_n2_dlubac;
    k_SDR_all{kk} = k_avg_best;
    k_ind_all{kk} = avgK_n2_ind;
    zMat_all{kk} = vertcat(zMat{:});
    KMat_all{kk} = vertcat(KMat{:});
    

    
end

kMaurer = vertcat(k_Maurer_all{:});
kDlubac = vertcat(k_Dlubac_all{:});
kInd = vertcat(k_ind_all{:});
kSDR = vertcat(k_SDR_all{:});
zMatTot = vertcat(zMat_all{:});
KMatTot = vertcat(KMat_all{:});

k_est_best = [kMaurer kInd kDlubac];
k_avg_best = [kMaurer kInd kDlubac];

figure(20)
plotKestKdpp_wdepth(KMatTot,k_est_best,k_avg_best,zMatTot,k_names)   
set(gca,'FontSize',16)

k_est_best = [kDlubac,kSDR];
k_avg_best = [kDlubac,kSDR];

figure(21)
plotKestKdpp_wdepth(KMatTot,k_est_best,k_avg_best,zMatTot,k_names)   
set(gca,'FontSize',16)

%legend(k_names_n2);

%save('bestKModels.mat','dlubacModel','bestSDRModel')
