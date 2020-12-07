clear

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load SDR_bestFit_1101_m1_n2.mat

sites = siteList;
SDR_b = squeeze(totalbMatrix(1,1,:))';
SDR_n = n;
SDR_n = repmat(n,length(siteList),1);
SDR_m = squeeze(totalmMatrix(1,1,:))';

%load Seevers_bestFit_1101_m1_n2_T2BAvg.mat
load Seevers_bestFit_0422_m1_n2_T2BMap.mat
sites = siteList;
Seevers_b = squeeze(totalbMatrix(1,1,:))';
Seevers_n = n;
Seevers_n = repmat(n,length(siteList),1);
Seevers_m = squeeze(totalmMatrix(1,1,:))';

% KGM (has been checked, looks good)
% KGM_tau = [1 1.7378 2.7227 3.4674];
% KGM_rho = [2.0606e-05 5.5335e-05 0.0077 1];
% KGM_m = [0 0 0 0];

% % New KGM with good T2B and m = 1 11/6/19
% KGM_tau = [1 1 1.135 1.5668];
% KGM_rho = [4.3251e-05 6.7421e-05 5.6754e-04 1.00e+02];
% KGM_m = [1 1 1 1];

% New KGM with map T2B=5.551 deg C and m = 1 4/8/20
KGM_tau = [1.0116 1 1.0965 1.4791];
KGM_rho = [7.8343e-05 4.7753e-05 6.9183e-04 100];
KGM_m = [1 1 1 1];

% SOE
SOE_n = [1 1 1 1];
SOE_b = [0.0045 0.0052 0.0158 0.0092];
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% K estimates
SDR_K = @(b,m,n,phi,T2ML) b.*(phi.^m).*(T2ML).^n;
SOE_K = @(b,n,SOE) b.*(SOE).^n;
Seevers_K = @(b,m,n,T2ML,T2B,phi) b.*(phi).^m.*((T2ML.^(-1) - T2B.^(-1)).^(-1)).^n;

Temp = 5.551;  % temperature in degress C 
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

for kk = 1:length(sites)
    [T2dist, T2logbins, nmrName] = loadRawNMRdata(sites{kk});
   
    [d, K, T2ML, phi, z, SumEch, logK, logT2ML, logPhi, SumEch_3s, SumEch_twm, ...
    SumEch_twm_3s] = loadnmrdata2(nmrName);

    SDR_Kest{kk} = SDR_K(SDR_b(kk),SDR_m(kk),SDR_n(kk),phi,T2ML);
    SOE_Kest{kk} = SOE_K(SOE_b(kk),SOE_n(kk),SumEch);
    
    Seevers_Kest{kk} = Seevers_K(Seevers_b(kk),Seevers_m(kk),Seevers_n(kk),...
        T2ML,2.2,phi);
    
    KGM_lKest{kk} = KGM_lK(KGM_rho(kk),KGM_tau(kk),KGM_m(kk),log10(phi),T2ML);
    KGM_Kest{kk} = 10.^KGM_lKest{kk};
    
    Ktot{kk} = K;
       
end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load SDR_maurer_bestFit_1101_m1_n2.mat
sites = siteList;
SDR_b = squeeze(totalbMatrix(1,1,:))';
SDR_n = n;
SDR_n = repmat(n,length(siteList),1);
SDR_m = squeeze(totalmMatrix(1,1,:))';

load Seevers_maurer_bestFit_0422_m1_n2_T2BMap.mat
sites = siteList;
Seevers_b = squeeze(totalbMatrix(1,1,:))';
Seevers_n = n;
Seevers_n = repmat(n,length(siteList),1);
Seevers_m = squeeze(totalmMatrix(1,1,:))';
% 
% % KGM (has been checked, looks good)
% KGM_tau = [1.6788 1.9498 1 1 1.1482 1.5311 1.4289 1.4125 1 1.3335];
% KGM_rho = [0.0100 4.2170e-4 5.2723e-05 9.5499e-05 1.4191e-04 100 100 100 1.6141e-04 1.4191e-04];
% KGM_m = [0 0 0 0 0 0 0 0 0 0];

%New KGM with T2B for T = 12.5 C and m = 1 11/6/19
% KGM_tau = [1 1.0471 1 1 1 1 1 1 1 1];
% KGM_rho = [100 100 1.2853e-04 1.300e-03 3.4000e-03 1.00e+02 1.00e+02 1.00e+02 6.9183e-04 2.2570e-04];
% KGM_m = [1 1 1 1 1 1 1 1 1 1];

% New KGM with T2B for T = 11.1 C from USGS Map and m = 1 4/7/20
KGM_tau = [1 1 1 1 1 1 1 1 1 1];
KGM_rho = [100 100 1.5668e-04 100 100 100 100 100 100 4.6559e-04];
KGM_m = [1 1 1 1 1 1 1 1 1 1];

% SOE
SOE_n = [1 1 1 1 1 1 1 1 1 1];
SOE_b = [0.0052 0.0025 0.0010 0.0033 0.0032 0.0080 0.0077 0.0074 0.0033 0.0013];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% K estimates
SDR_K = @(b,m,n,phi,T2ML) b.*(phi.^m).*(T2ML).^n;
SOE_K = @(b,n,SOE) b.*(SOE).^n;
Seevers_K = @(b,m,n,T2ML,T2B,phi) b.*(phi).^m.*((T2ML.^(-1) - T2B.^(-1)).^(-1)).^n;

Temp = 5.551;  % temperature in degress C 
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

for kk = 1:length(sites)
    [T2dist, T2logbins, nmrName] = loadRawNMRdata(sites{kk});
   
    [d, K, T2ML, phi, z, SumEch, logK, logT2ML, logPhi, SumEch_3s, SumEch_twm, ...
    SumEch_twm_3s] = loadnmrdata2(nmrName);

    SDR_Kest_Maurer{kk} = SDR_K(SDR_b(kk),SDR_m(kk),SDR_n(kk),phi,T2ML);
    SOE_Kest_Maurer{kk} = SOE_K(SOE_b(kk),SOE_n(kk),SumEch);
    
    Seevers_Kest_Maurer{kk} = Seevers_K(Seevers_b(kk),Seevers_m(kk),Seevers_n(kk),...
        T2ML,2.2,phi);
    
    KGM_lKest_Maurer{kk} = KGM_lK(KGM_rho(kk),KGM_tau(kk),KGM_m(kk),log10(phi),T2ML);
    KGM_Kest_Maurer{kk} = 10.^KGM_lKest_Maurer{kk};
    
    Ktot_Maurer{kk} = K;
       
end

K_all = vertcat(Ktot{:}, Ktot_Maurer{:});
k_estimates = [vertcat(SDR_Kest{:}, SDR_Kest_Maurer{:}) vertcat(SOE_Kest{:},SOE_Kest_Maurer{:})...
    vertcat(Seevers_Kest{:},Seevers_Kest_Maurer{:}) vertcat(KGM_Kest{:},KGM_Kest_Maurer{:})];
k_names = {'SDR','SOE','Seevers','KGM'};

[allSigns kDiffFactorTemp] = estimateKdiffFactor_withSign(K_all, k_estimates, 1);
SDR_Kdiff = log10(kDiffFactorTemp(1:end,1)).*allSigns(1:end,1);
SOE_Kdiff = log10(kDiffFactorTemp(1:end,2)).*allSigns(1:end,2);
Seevers_Kdiff = log10(kDiffFactorTemp(1:end,3)).*allSigns(1:end,3);
KGM_Kdiff = log10(kDiffFactorTemp(1:end,4)).*allSigns(1:end,4);


Khigh = K_all(K_all >5*10^(-4));
KDiffFactorSDRHigh = SDR_Kdiff(K_all > 5*10^(-4));
KDiffFactorSOEHigh = SOE_Kdiff(K_all > 5*10^(-4));
KDiffFactorSeeversHigh = Seevers_Kdiff(K_all > 5*10^(-4));
KDiffFactorKGMHigh = KGM_Kdiff(K_all > 5*10^(-4));
figure(1)
KmodelDiffHist(KDiffFactorSDRHigh, KDiffFactorSOEHigh, KDiffFactorSeeversHigh, KDiffFactorKGMHigh)
suptitle('"High"-K (K > 5*10^-4 m/s)')

Kmed = K_all(K_all > 10^(-5) & K_all < 5*10^(-4));
KDiffFactorSDRMed = SDR_Kdiff(K_all > 10^(-5) & K_all < 5*10^(-4));
KDiffFactorSOEMed = SOE_Kdiff(K_all > 10^(-5) & K_all < 5*10^(-4));
KDiffFactorSeeversMed = Seevers_Kdiff(K_all > 10^(-5) & K_all < 5*10^(-4));
KDiffFactorKGMMed = KGM_Kdiff(K_all > 10^(-5) & K_all < 5*10^(-4));
figure(2)
KmodelDiffHist(KDiffFactorSDRMed, KDiffFactorSOEMed, KDiffFactorSeeversMed, KDiffFactorKGMMed)
suptitle('"Medium"-K (K > 10^-5 m/s & K < 5*10^-4 m/s)')

Klow = K_all(K_all < 10^(-5));
KDiffFactorSDRLow = SDR_Kdiff(K_all < 10^(-5));
KDiffFactorSOELow = SOE_Kdiff(K_all < 10^(-5));
KDiffFactorSeeversLow = Seevers_Kdiff(K_all < 10^(-5));
KDiffFactorKGMLow = KGM_Kdiff(K_all < 10^(-5));
figure(3)
KmodelDiffHist(KDiffFactorSDRLow, KDiffFactorSOELow, KDiffFactorSeeversLow, KDiffFactorKGMLow)
suptitle('"Low"-K (K < 10^-5 m/s)')


