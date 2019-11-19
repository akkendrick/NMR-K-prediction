clear
close all

sites = [{'Site1-WellG6'} {'Site1-WellG5'} {'Site2-WellPN1'} {'Site2-WellPN2'}];

% Best fit params from 6/18/19
% SDR
SDR_b = [0.0024 0.0040 0.0075 0.0047];
SDR_m = [0 0 0 0];
SDR_n = [2 2 2 2];

% Seevers
Seevers_b = [0.0016 0.0032 0.0058 0.0038];
Seevers_n = [2 2 2 2];
Seevers_m = [0 0 0 0];

% KGM (has been checked, looks good)
KGM_tau = [1 1.7378 2.7227 3.4674];
KGM_rho = [2.0606e-05 5.5335e-05 0.0077 1];
% KGM_tau = [1 1.7378 2.7227 2.5];
% KGM_rho = [2.0606e-05 5.5335e-05 10 1]; 

KGM_m = [0 0 0 0];

% SOE
SOE_n = [1 1 1 1];
SOE_b = [0.0045 0.0052 0.0158 0.0092];

% K estimates
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

for kk = 1:length(sites)
    [T2dist, T2logbins, nmrName] = loadRawNMRdata(sites{kk});
   
    [d, K, T2ML, phi, z, SumEch, logK, logT2ML, logPhi, SumEch_3s, SumEch_twm, ...
    SumEch_twm_3s] = loadnmrdata2(nmrName);

    SDR_Kest = SDR_K(SDR_b(kk),SDR_m(kk),SDR_n(kk),phi,T2ML);
    SOE_Kest = SOE_K(SOE_b(kk),SOE_n(kk),SumEch);
    
    Seevers_Kest = Seevers_K(Seevers_b(kk),Seevers_m(kk),Seevers_n(kk),...
        T2ML,2.2,phi);
    
    KGM_lKest = KGM_lK(KGM_rho(kk),KGM_tau(kk),KGM_m(kk),log10(phi),T2ML);
    KGM_Kest = 10.^KGM_lKest;
    
    k_estimates = [SDR_Kest SOE_Kest Seevers_Kest KGM_Kest];
    k_names = {'SDR','SOE','Seevers','KGM'};
    k_sym = {'+','o','*','x','^','v'};
    
    kDiffFactorTemp = estimateKdiffFactor(K, k_estimates, 1);
    kDiffFactorMean{:,kk} = mean(kDiffFactorTemp);
    kDiffFactor{:,kk} = kDiffFactorTemp;
    
    figure(1)
    subplot(4,1,kk)
    
    plotKestKdpp_v2(K,k_estimates,k_estimates,k_names)
    title(sites{kk})
    
end
    
    