% Plot basic K, T2ml, water table, phi profiles

clear
close all

load 'SDR_bestFit_table.mat'
load 'Seevers_bestFit_table.mat'
load 'TC_bestFit_240_table.mat'

sites = {'Site1-WellG5sep','Site1-WellG6sep','Site2-WellPN1','Site2-WellPN2'};

waterTable = [2.0469,2.1248,5.0285,4.7476]; % rel ground surface
%waterTable = [3.004,2.963,5.727,5.408]; % rel top of casing
depthOffsets = [0.75,0.95,0.75,0.75];


% Sonic depths are relative to ground surface
load('sonicCoreT2B.mat','sonicCoreT2BData')
T2B_depth = sonicCoreT2BData.Depthm;
T2B_depth = flipud(T2B_depth);
 
T2B_peak = sonicCoreT2BData.T2Bpeak;
T2B_peak = flipud(T2B_peak); 

model = 'bootstrap';

%% Define all K models
SDR_K = @(b,m,n,phi,T2ML) b.*(phi.^m).*(T2ML).^n;
SOE_K = @(b,n,SOE) b.*(SOE).^n;
Seevers_K = @(b,m,n,T2ML,T2B,phi) b.*(phi).^m.*((T2ML.^(-1) - T2B.^(-1)).^(-1)).^n;

lkTC = @(c,m,n,lPhi,logFrac) log10(c) + m*lPhi + n*(logFrac);


Temp = 20;  % temperature in degress C 
rho = @(Tt) 1000*(1 - ((Tt+288.94)./(508929*(Tt+68.12))).*(Tt-3.98).^2); % kg/m^3
eta = @(Tt) 0.0013 - 1.7e-5*Tt;         % Pa -s
%Tb = @(Tt) 3.3 + 0.044*(Tt - 35);       % seconds

Tb= @(Tt) 3.3 + 0.044.*(Tt - 35);                  % Bulk relaxation time constant (s)
D = @(Tt) (1.0413 + 0.039828*Tt + 0.00040318*Tt.^2).*1e-9;  % m^2/s 
g = 9.8;    %m/s^2
tort = 1/(1.5^2); 
t1 = (rho(Temp)*g)/(8*eta(Temp)); % 

num2 = @(T2) 4*D(Temp)*Tb(Temp)*T2;
denom2 = @(T2) Tb(Temp) - T2; 
    
f12 = @(rho) (D(Temp)./rho);  
SQterm = @(rho,T2) sqrt(f12(rho).^2 + (num2(T2)./denom2(T2))); 

KGM_lK = @(rho,tau,m,lphi,T2) log10(1/tau^2) + log10(t1) + m*lphi + 2*log10(SQterm(rho,T2)-f12(rho)); 
lK_kgm = @(A, B, T, T2, phi) KGMfun([A, B], [T2]); 

%%

for kk = 1:length(sites)
    site = sites{kk};
    baseDir = '/Volumes/GoogleDrive/My Drive/USGS Project/USGS Data/';
    baseNMRDir = '/Volumes/GoogleDrive/My Drive/USGS Project/NMR-K-prediction/';
 
 %    baseDir = 'I:\My Drive\USGS Project\USGS Data\';
%    baseNMRDir = 'I:\My Drive\USGS Project\NMR-K-prediction\';
    
    T2dist = [];
    T2logbins = [];
    dparam = [];
    kk
    clear kTC lK KMat T2MLMat phiMat zMat SumEchMat logKMat logT2MLMat logPhiMat NMRphi kSDR interpT2B kSOE kKGM
    
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

        indexQuotient{:,1} = TC_indexQuotient{4}; 
        
        nmrDepths = z(1:end,1);
        interpT2B = interp1(T2B_depth, T2B_peak, nmrDepths)*10^-3; %Convert from ms to s
        
        siteListInd = 4;
        if strcmp(model,'boostrap')
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
        
        indexQuotient{:,1} = TC_indexQuotient{5}; 
        indexQuotient{:,2} = TC_indexQuotient{6}; 

        
        siteListInd = [5 6];
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

        indexQuotient{:,1} = TC_indexQuotient{1}; 

        
        siteListInd = [1];
        if strcmp(model,'bootstrap')
            SDRb = SDR_Bootstrap(1,siteListInd);
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

        indexQuotient{:,1} = TC_indexQuotient{2}; 
        indexQuotient{:,2} = TC_indexQuotient{3}; 

        siteListInd = [2 3];
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
        
        indexQuotient{:,1} = TC_indexQuotient{7}; 


        siteListInd = [7];
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
        indexQuotient{:,1} = TC_indexQuotient{8}; 
        siteListInd = [8];
        
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
    end
    
    %% Make a matrix of all models so far
    
    noPorosity = 1;
    
    for j = 1:length(SDRn)
        kSDR{j} = SDR_K(SDRb(j),SDRm(j),SDRn(j),phiMat{j},T2MLMat{j});
        kSOE{j} = SOE_K(SOEb(j),SOEn(j),SumEchMat{j});
        
        lkTC_est{j} = lkTC(TC33b(j),TC33m(j),TC33n(j),logPhiMat{j},indexQuotient{j});
        kTC{j} = 10.^lkTC_est{j};
        
        A(j) = 2*log10(1/KGMtau(j));
        
        lK{j} = lK_kgm(A(j),KGMrho(j),20,T2MLMat{j},logPhiMat{j});
        kKGM{j} = 10.^lK{j};
%         lK_kgm = @(A, B, T, T2, phi) KGMfun([A, B], [T2]); 

%         kKGM{j} = KGM_lK(KGMrho(j),KGMtau(j),2,logPhiMat{j},T2MLMat{j});
%         kKGM{j} = 10.^kKGM{j};
        
        if ~isempty(interpT2B)
            kSeevers{:,j} = Seevers_K(Seeversb(j),Seeversm(j),Seeversn(j),T2MLMat{j},interpT2B{j},phiMat{j});
            %kSDR_max{:,j} = SDR_K(SDRb(j),SDRm(j),SDRn(j),phiMat{j},interpT2B{j});
            kSDR_max = {};
            %k_names = [{'DPP'} {'SDR'} {'SOE'} {'KGM'} {'TC 240'} {'Seevers'} {'SDR Max'}];
            k_names = [{'DPP'} {'SDR'} {'SOE'} {'KGM'} {'TC 240'} {'Seevers'} {'Water Table'}];

        else
            kSeevers = {};
            kSDR_max = {};

            k_names = [{'DPP'} {'SDR'} {'SOE'} {'KGM'} {'TC 240'} {'Water Table'}];

        end
        
    %k_estimates = [vertcat(k_bootstrap{:}) vertcat(k_direct{:}) vertcat(k_mcmc{:})];
    %k_names = [{'DPP'} {'Bootstrap'} {'Direct'} {'MCMC'}];
    %k_sym = [{'+'} {'+'} {'+'}]
    %T2dist = flip(T2dist);


    end
    
    if ~isempty(interpT2B)
        errorMatrix(1,kk) = computeError(vertcat(KMat{:}),vertcat(kSDR{:}));
        errorMatrix(2,kk) = computeError(vertcat(KMat{:}),vertcat(kSOE{:}));
        errorMatrix(3,kk) = computeError(vertcat(KMat{:}),vertcat(kKGM{:}));
        errorMatrix(4,kk) = computeError(vertcat(KMat{:}),vertcat(kTC{:}));
        errorMatrix(5,kk) = computeError(vertcat(KMat{:}),vertcat(kSeevers{:}));
        
    else
        errorMatrix(1,kk) = computeError(vertcat(KMat{:}),vertcat(kSDR{:}));
        errorMatrix(2,kk) = computeError(vertcat(KMat{:}),vertcat(kSOE{:}));
        errorMatrix(3,kk) = computeError(vertcat(KMat{:}),vertcat(kKGM{:}));
        errorMatrix(4,kk) = computeError(vertcat(KMat{:}),vertcat(kTC{:}));
        
    end

    
    k_estimates = [vertcat(kSDR{:}) vertcat(kSOE{:}) vertcat(kKGM{:}) vertcat(kTC{:}) vertcat(kSeevers{:}) vertcat(kSDR_max{:})];
    k_sym = [{'+'} {'o'} {'x'} {'v'} {'>'},{'>'}];
    
    plotKPhiwithDepth(vertcat(KMat{:}),NMRphi,waterTable(kk),vertcat(zMat{:}),T2dist,...
        T2logbins,k_estimates,k_names,k_sym,noPorosity)
    title(site)
    
    outName = strcat(site,'_profile.fig');
    savefig(outName)
    
    outName = strcat(site,'_profile.png');
    print(outName,'-dpng')
    
    
    
    plotKestKdpp(vertcat(KMat{:}),k_estimates,k_names,k_sym)   
    title(site)  
    
    outName = strcat(site,'_Kcomp.fig');
    savefig(outName)
    
    outName = strcat(site,'_Kcomp.png');
    print(outName,'-dpng')
    
    kDiffFactor = estimateKdiffFactor(vertcat(KMat{:}),k_estimates);
    
    display('kDiffFactor')
    meanDiffFactor = mean(kDiffFactor)

end