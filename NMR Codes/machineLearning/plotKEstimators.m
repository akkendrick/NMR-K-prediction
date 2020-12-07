% Plot K estimators (T2ML, Sum of Echoes) to estimate when models are more
% appropriate

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sites = [{'Site1-WellG6'} {'Site1-WellG5'} {'Site2-WellPN1'} {'Site2-WellPN2'}];

% Best fit params from 6/18/19
% SDR
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

% New KGM with good T2B and m = 1 11/6/19
KGM_tau = [1 1 1.135 1.5668];
KGM_rho = [4.3251e-05 6.7421e-05 5.6754e-04 1.00e+02];
KGM_m = [1 1 1 1];

% SOE
SOE_n = [1 1 1 1];
SOE_b = [0.0045 0.0052 0.0158 0.0092];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sites_Maurer = {'dpnmr_larned_east','dpnmr_larned_lwph','dpnmr_larned_west',...
   'dpnmrA11','dpnmrA12','dpnmrC1S','dpnmrC1SE','dpnmrC1SW',...
   'dpnmr_leque_east','dpnmr_leque_west'}

load SDR_maurer_bestFit_1101_m1_n2.mat
sites_Maurer = siteList;
SDR_b_Maurer = squeeze(totalbMatrix(1,1,:))';
SDR_n_Maurer = n;
SDR_n_Maurer = repmat(n,length(siteList),1);
SDR_m_Maurer = squeeze(totalmMatrix(1,1,:))';

load Seevers_maurer_bestFit_1101_m1_n2_T2B2.31.mat
sites_Maurer = siteList;
Seevers_b_Maurer = squeeze(totalbMatrix(1,1,:))';
Seevers_n_Maurer = n;
Seevers_n_Maurer = repmat(n,length(siteList),1);
Seevers_m_Maurer = squeeze(totalmMatrix(1,1,:))';

%New KGM with T2B for T = 12.5 C and m = 1 11/6/19
KGM_tau_Maurer = [1 1.0471 1 1 1 1 1 1 1 1];
KGM_rho_Maurer = [100 100 1.2853e-04 1.300e-03 3.4000e-03 1.00e+02 1.00e+02 1.00e+02 6.9183e-04 2.2570e-04];
KGM_m = [1 1 1 1 1 1 1 1 1 1];

% SOE
SOE_n = [1 1 1 1 1 1 1 1 1 1];
SOE_b_Maurer = [0.0052 0.0025 0.0010 0.0033 0.0032 0.0080 0.0077 0.0074 0.0033 0.0013];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sonic depths are relative to ground surface
load('sonicCoreT2B.mat','sonicCoreT2BData')
T2B_depth = sonicCoreT2BData.Depthm;
T2B_depth = flipud(T2B_depth);
 
T2B_peak = sonicCoreT2BData.T2Bpeak;
T2B_peak = flipud(T2B_peak); 
goodT2B = T2B_peak(T2B_peak > 1000);
T2Bavg = mean(goodT2B)*10^-3;

SDR_bCorrFactor = mean(SDR_b_Maurer)/mean(SDR_b);
SOE_bCorrFactor = mean(SOE_b_Maurer)/mean(SOE_b);
Seevers_bCorrFactor = mean(Seevers_b_Maurer)/mean(Seevers_b);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for kk = 1:length(sites)
    [T2dist, T2logbins, nmrName] = loadRawNMRdata(sites{kk});
   
    [d, K, T2ML, phi, z, SumEch, logK, logT2ML, logPhi, SumEch_3s, SumEch_twm, ...
    SumEch_twm_3s] = loadnmrdata2(nmrName);
    
    Kall{kk} = K;
    Phiall{kk} = phi;
    T2MLall{kk} = T2ML;
    SumEchAll{kk} = SumEch;
   
end

for kk = 1:length(sites_Maurer)
    [T2dist, T2logbins, nmrName] = loadRawNMRdata(sites_Maurer{kk});
   
    [d, K_Maurer, T2ML_Maurer, phi_Maurer, z, SumEch_Maurer, logK, logT2ML, logPhi, SumEch_3s, SumEch_twm, ...
    SumEch_twm_3s] = loadnmrdata2(nmrName);
    
    Kall_Maurer{kk} = K_Maurer;
    Phiall_Maurer{kk} = phi_Maurer;
    T2MLall_Maurer{kk} = T2ML_Maurer;
    SumEchAll_Maurer{kk} = SumEch_Maurer;
        
end

Kall_vec = vertcat(Kall{:});
Kall_Maurer_vec = vertcat(Kall_Maurer{:});
T2MLall_vec = vertcat(T2MLall{:});
T2MLall_Maurer_vec = vertcat(T2MLall_Maurer{:});
SumEchAll_vec = vertcat(SumEchAll{:});
SumEchAll_Maurer_vec = vertcat(SumEchAll_Maurer{:});

%%
% Plot T2ML
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
box on
grid on
hold on

corrT2ML2all_vec = (T2MLall_vec.^2)/SDR_bCorrFactor;

%scatter(T2MLall_vec.^2, Kall_vec, 40,'Filled')
scatter(corrT2ML2all_vec, Kall_vec, 40,'Filled')
scatter(T2MLall_Maurer_vec.^2,Kall_Maurer_vec,40,'Filled')

set(gca,'XScale','log')
set(gca,'YScale','log')
set(gca,'FontSize',16)

ylim([10^-7,10^-2])
xlim([10^-4, 10^0])

ylabel('Hydraulic Conductivity (m/s)')
xlabel('T2ML.^2 (s)')

legend({'Wisc','Maurer'},'Location','northwest')

% Plot SOE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(2)
box on
grid on
hold on

corrSumEchAll_vec = SumEchAll_vec/SOE_bCorrFactor;

%scatter(SumEchAll_vec, Kall_vec, 40,'Filled')
scatter(corrSumEchAll_vec, Kall_vec, 40, 'Filled')
scatter(SumEchAll_Maurer_vec,Kall_Maurer_vec,40,'Filled')

set(gca,'XScale','log')
set(gca,'YScale','log')
set(gca,'FontSize',16)

ylim([10^-7,10^-2])

ylabel('Hydraulic Conductivity (m/s)')
xlabel('Sum of Spin-Echoes')

legend({'Wisc','Maurer'},'Location','northwest')

% Plot T2Seevers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T2B = 2.2293; % Avg pore water from Wisconsin
T2Seevers_vec = ((T2MLall_vec.^(-1) - T2B.^(-1)).^(-1));

corrT2Seevers2all_vec = (T2Seevers_vec.^2)/Seevers_bCorrFactor;

T2B = 2.3100; %s T2B from Dlugosch for 12.5 deg C
T2Seevers_Maurer_vec = ((T2MLall_Maurer_vec.^(-1) - T2B.^(-1)).^(-1));

figure(3)
box on
grid on
hold on

scatter(corrT2Seevers2all_vec, Kall_vec, 40, 'Filled')
scatter(T2Seevers_Maurer_vec.^2,Kall_Maurer_vec,40,'Filled')

set(gca,'XScale','log')
set(gca,'YScale','log')
set(gca,'FontSize',16)

ylim([10^-7,10^-2])
xlim([10^-4, 10^0])

ylabel('Hydraulic Conductivity (m/s)')
xlabel('T2Seevers^2 (s)')

legend({'Wisc','Maurer'},'Location','northwest')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(4)
subplot(2,1,1)
box on
grid on
hold on

%scatter(corrT2ML2all_vec, Kall_vec, 40,'Filled')
%scatter(corrT2Seevers2all_vec, Kall_vec, 40, 'Filled')

scatter(T2MLall_vec.^2, Kall_vec, 40,'Filled')
scatter(T2Seevers_vec.^2, Kall_vec, 40,'Filled')


set(gca,'XScale','log')
set(gca,'YScale','log')
set(gca,'FontSize',16)


subplot(2,1,2)
box on
grid on
hold on

scatter(T2MLall_Maurer_vec.^2, Kall_Maurer_vec, 40,'Filled')
scatter(T2Seevers_Maurer_vec.^2, Kall_Maurer_vec, 40, 'Filled')

set(gca,'XScale','log')
set(gca,'YScale','log')
set(gca,'FontSize',16)



