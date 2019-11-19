% Plot K diff Histogram
sites = [{'Site1-WellG6'} {'Site1-WellG5'} {'Site2-WellPN1'} {'Site2-WellPN2'}];

sites_Maurer = {'dpnmr_larned_east','dpnmr_larned_lwph','dpnmr_larned_west',...
   'dpnmrA11','dpnmrA12','dpnmrC1S','dpnmrC1SE','dpnmrC1SW',...
   'dpnmr_leque_east','dpnmr_leque_west'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SDR_K = @(b,m,n,phi,T2ML) b.*(phi.^m).*(T2ML).^n;
SOE_K = @(b,n,SOE) b.*(SOE).^n;
Seevers_K = @(b,m,n,T2ML,T2B,phi) b.*(phi).^m.*((T2ML.^(-1) - T2B.^(-1)).^(-1)).^n;

Temp = 7.4;  % temperature in degress C 
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

T2B = 2.2293;

load SDR_bestFit_1101_m1_n2.mat
sites = siteList;
SDR_b = squeeze(totalbMatrix(1,1,:))';
SDR_n = n;
SDR_n = repmat(n,length(siteList),1);
SDR_m = squeeze(totalmMatrix(1,1,:))';

load Seevers_bestFit_1101_m1_n2_T2BAvg.mat
sites = siteList;
Seevers_b = squeeze(totalbMatrix(1,1,:))';
Seevers_n = n;
Seevers_n = repmat(n,length(siteList),1);
Seevers_m = squeeze(totalmMatrix(1,1,:))';

% KGM (has been checked, looks good)
KGM_tau = [1 1 1.135 1.5668];
KGM_rho = [4.3251e-05 6.7421e-05 5.6754e-04 1.00e+02];
KGM_m = [1 1 1 1];

% SOE
SOE_n = [1 1 1 1];
SOE_b = [0.0045 0.0052 0.0158 0.0092];

for kk = 1:length(sites)

    
    [T2dist, T2logbins, nmrName] = loadRawNMRdata(sites{kk});
   
    T2dist_matrix{kk} = T2dist;
    T2logbins_matrix{kk} = T2logbins;
    
    [d, K, T2ML, phi, z, SumEch, logK, logT2ML, logPhi, SumEch_3s, SumEch_twm, ...
    SumEch_twm_3s] = loadnmrdata2(nmrName);
    
    Kall{kk} = K;
    Phiall{kk} = phi;
    T2MLall{kk} = T2ML;
    SumEchAll{kk} = SumEch;
    
    SDR_Kest{kk} = SDR_K(SDR_b(kk),SDR_m(kk),SDR_n(kk),phi,T2ML);
    SOE_Kest{kk} = SOE_K(SOE_b(kk),SOE_n(kk),SumEch);
    
    Seevers_Kest{kk} = Seevers_K(Seevers_b(kk),Seevers_m(kk),Seevers_n(kk),...
        T2ML,T2B,phi);
    
    KGM_lKest{kk} = KGM_lK(KGM_rho(kk),KGM_tau(kk),KGM_m(kk),log10(phi),T2ML);
    KGM_Kest{kk} = 10.^KGM_lKest{kk};
    
    Ktot{kk} = K;
    
    [SDR_sign{kk} SDR_diffFactor{kk}] = estimateKdiffFactor_withSign(K,SDR_Kest{kk},1);
    [Seevers_sign{kk} Seevers_diffFactor{kk}] = estimateKdiffFactor_withSign(K,Seevers_Kest{kk},1);
    [SOE_sign{kk} SOE_diffFactor{kk}] = estimateKdiffFactor_withSign(K,SOE_Kest{kk},1);
    [KGM_sign{kk} KGM_diffFactor{kk}] = estimateKdiffFactor_withSign(K,KGM_Kest{kk},1);
    
end

T2B = 2.31;

load SDR_maurer_bestFit_1101_m1_n2.mat

sites = siteList;
SDR_b = squeeze(totalbMatrix(1,1,:))';
SDR_n = n;
SDR_n = repmat(n,length(siteList),1);
SDR_m = squeeze(totalmMatrix(1,1,:))';

load Seevers_maurer_bestFit_1101_m1_n2_T2B2.31.mat
sites = siteList;
Seevers_b = squeeze(totalbMatrix(1,1,:))';
Seevers_n = n;
Seevers_n = repmat(n,length(siteList),1);
Seevers_m = squeeze(totalmMatrix(1,1,:))';

% KGM (has been checked, looks good)
KGM_tau = [1 1.0471 1 1 1 1 1 1 1 1];
KGM_rho = [100 100 1.2853e-04 1.300e-03 3.4000e-03 1.00e+02 1.00e+02 1.00e+02 6.9183e-04 2.2570e-04];
KGM_m = [1 1 1 1 1 1 1 1 1 1];

% SOE
SOE_n = [1 1 1 1 1 1 1 1 1 1];
SOE_b = [0.0052 0.0025 0.0010 0.0033 0.0032 0.0080 0.0077 0.0074 0.0033 0.0013];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Temp = 12.5;  % temperature in degress C 
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for kk = 1:length(sites_Maurer)
    [T2dist, T2logbins, nmrName] = loadRawNMRdata(sites_Maurer{kk});
   
    T2dist_matrix_Maurer{kk} = T2dist;
    T2logbins_matrix_Maurer{kk} = T2logbins;
    
    [d, K_Maurer, T2ML_Maurer, phi_Maurer, z, SumEch_Maurer, logK, logT2ML, logPhi, SumEch_3s, SumEch_twm, ...
    SumEch_twm_3s] = loadnmrdata2(nmrName);
    
    Kall_Maurer{kk} = K_Maurer;
    Phiall_Maurer{kk} = phi_Maurer;
    T2MLall_Maurer{kk} = T2ML_Maurer;
    
    SumEchAll_Maurer{kk} = SumEch_Maurer;
    
    SDR_Kest_Maurer{kk} = SDR_K(SDR_b(kk),SDR_m(kk),SDR_n(kk),phi_Maurer,T2ML_Maurer);
    SOE_Kest_Maurer{kk} = SOE_K(SOE_b(kk),SOE_n(kk),SumEch_Maurer);

    Seevers_Kest_Maurer{kk} = Seevers_K(Seevers_b(kk),Seevers_m(kk),Seevers_n(kk),...
    T2ML_Maurer,T2B,phi_Maurer);

    KGM_lKest{kk} = KGM_lK(KGM_rho(kk),KGM_tau(kk),KGM_m(kk),log10(phi_Maurer),T2ML_Maurer);
    KGM_Kest_Maurer{kk} = 10.^KGM_lKest{kk};

    Ktot_Maurer{kk} = K_Maurer;

    [SDR_sign_Maurer{kk} SDR_diffFactor_Maurer{kk}] = estimateKdiffFactor_withSign(K_Maurer,SDR_Kest_Maurer{kk},1);
    [Seevers_sign_Maurer{kk} Seevers_diffFactor_Maurer{kk}] = estimateKdiffFactor_withSign(K_Maurer,Seevers_Kest_Maurer{kk},1);
    [SOE_sign_Maurer{kk} SOE_diffFactor_Maurer{kk}] = estimateKdiffFactor_withSign(K_Maurer,SOE_Kest_Maurer{kk},1);
    [KGM_sign_Maurer{kk} KGM_diffFactor_Maurer{kk}] = estimateKdiffFactor_withSign(K_Maurer,KGM_Kest_Maurer{kk},1);
        
end

Kall_vec = vertcat(Kall{:});
Kall_Maurer_vec = vertcat(Kall_Maurer{:});
T2MLall_vec = vertcat(T2MLall{:});
T2MLall_Maurer_vec = vertcat(T2MLall_Maurer{:});

SDR_diffFactor_all = vertcat(SDR_diffFactor{:});
Seevers_diffFactor_all = vertcat(Seevers_diffFactor{:});
SOE_diffFactor_all = vertcat(SOE_diffFactor{:});
KGM_diffFactor_all = vertcat(KGM_diffFactor{:});

SDR_diffFactor_all_Maurer = vertcat(SDR_diffFactor_Maurer{:});
Seevers_diffFactor_all_Maurer = vertcat(Seevers_diffFactor_Maurer{:});
SOE_diffFactor_all_Maurer = vertcat(SOE_diffFactor_Maurer{:});
KGM_diffFactor_all_Maurer = vertcat(KGM_diffFactor_Maurer{:});

SDR_sign_all = vertcat(SDR_sign{:});
Seevers_sign_all = vertcat(Seevers_sign{:});
SOE_sign_all = vertcat(SOE_sign{:});
KGM_sign_all = vertcat(KGM_sign{:});

SDR_sign_all_Maurer = vertcat(SDR_sign_Maurer{:});
Seevers_sign_all_Maurer = vertcat(Seevers_sign_Maurer{:});
SOE_sign_all_Maurer = vertcat(SOE_sign_Maurer{:});
KGM_sign_all_Maurer = vertcat(KGM_sign_Maurer{:});

SDR_diffFactor_total = [SDR_diffFactor_all; SDR_diffFactor_all_Maurer];
Seevers_diffFactor_total = [Seevers_diffFactor_all; Seevers_diffFactor_all_Maurer];
SOE_diffFactor_total = [SOE_diffFactor_all; SOE_diffFactor_all_Maurer];
KGM_diffFactor_total = [KGM_diffFactor_all; KGM_diffFactor_all_Maurer];

SDR_sign_total = [SDR_sign_all; SDR_sign_all_Maurer];
Seevers_sign_total = [Seevers_sign_all; Seevers_sign_all_Maurer];
SOE_sign_total = [SOE_sign_all; SOE_sign_all_Maurer];
KGM_sign_total = [KGM_sign_all; KGM_sign_all_Maurer];

% SDR_diffFactor_total = [SDR_diffFactor_all];
% Seevers_diffFactor_total = [Seevers_diffFactor_all];
% SOE_diffFactor_total = [SOE_diffFactor_all];
% KGM_diffFactor_total = [KGM_diffFactor_all];
% 
% SDR_sign_total = [SDR_sign_all];
% Seevers_sign_total = [Seevers_sign_all];
% SOE_sign_total = [SOE_sign_all];
% KGM_sign_total = [KGM_sign_all];


% Manually assign diff factor for a factor above 20
SDR_diffFactor_total(SDR_diffFactor_total >= 100) = 18;
Seevers_diffFactor_total(Seevers_diffFactor_total >= 100) = 18;
SOE_diffFactor_total(SOE_diffFactor_total >= 100) = 18;
KGM_diffFactor_total(KGM_diffFactor_total >= 100) = 18;

%% 
% Plot histogram of Kdiff Factor
%edges = [1 1.5 2 2.5 3 4 5 10 20 40 80];
edgesLog = [-2 -1 -0.8 -0.6 -0.4 -0.2 0 0.2 0.4 0.6 0.8 1 2];
edges = 10.^edgesLog;

figure(4)

subplot(2,2,1)
hold on
grid on
box on
grid minor

title('SDR')

histDataSDR = SDR_sign_total.*log10(SDR_diffFactor_total);

%histogram(SDR_diffFactor_all, edges)
histogram(histDataSDR,edgesLog)
ylim([0 50])
%ylim([0 30])

xticklabels({'-100','-10','0','10','100'})
xlabel('K Difference Factor')
ylabel('Counts')
set(gca,'FontSize',14)

subplot(2,2,2)
hold on
grid on
box on
grid minor

title('Seevers')
histDataSeevers = Seevers_sign_total.*log10(Seevers_diffFactor_total);

%histogram(Seevers_diffFactor_all, edges)
histogram(histDataSeevers,edgesLog)
ylim([0 50])
%ylim([0 30])

xticklabels({'-100','-10','0','10','100'})
xlabel('K Difference Factor')
ylabel('Counts')
set(gca,'FontSize',14)


subplot(2,2,3)

hold on
grid on
box on
grid minor

title('SOE')
histDataSOE = SOE_sign_total.*log10(SOE_diffFactor_total);

%histogram(SOE_diffFactor_all, edges)
histogram(histDataSOE,edgesLog)
ylim([0 50])
%ylim([0 30])

xticklabels({'-100','-10','0','10','100'})
xlabel('K Difference Factor')
ylabel('Counts')
set(gca,'FontSize',14)


subplot(2,2,4)
hold on
grid on
box on
grid minor

title('KGM')
histDataKGM = KGM_sign_total.*log10(KGM_diffFactor_total);


%histogram(KGM_diffFactor_all, edges)
histogram(histDataKGM,edgesLog)
ylim([0 50])
%ylim([0 25])

xticklabels({'-100','-10','0','10','100'})
xlabel('K Difference Factor')
ylabel('Counts')
set(gca,'FontSize',14)

