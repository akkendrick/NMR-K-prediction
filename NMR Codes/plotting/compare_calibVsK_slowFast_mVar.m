% Compare fast/slow diffusion, different m coeff

% Estimate calibration as a function of K

clear
%close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
KGM_m = [0 0 0 0];

% SOE
SOE_n = [1 1 1 1];
SOE_b = [0.0045 0.0052 0.0158 0.0092];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sites_Maurer = {'dpnmr_larned_east','dpnmr_larned_lwph','dpnmr_larned_west',...
   'dpnmrA11','dpnmrA12','dpnmrC1S','dpnmrC1SE','dpnmrC1SW',...
   'dpnmr_leque_east','dpnmr_leque_west'}

load SDR_maurer_bestFit_1015_m0_n2.mat

sites_Maurer = siteList;
SDR_b_Maurer = squeeze(totalbMatrix(1,1,:))';
SDR_n_Maurer = n;
SDR_n_Maurer = repmat(n,length(siteList),1);
SDR_m_Maurer = squeeze(totalmMatrix(1,1,:))';

load Seevers_maurer_bestFit_1015_m0_n2.mat
sites_Maurer = siteList;
Seevers_b_Maurer = squeeze(totalbMatrix(1,1,:))';
Seevers_n_Maurer = n;
Seevers_n_Maurer = repmat(n,length(siteList),1);
Seevers_m_Maurer = squeeze(totalmMatrix(1,1,:))';

% KGM (has been checked, looks good)
KGM_tau_Maurer = [1.6788 1.9498 1 1 1.1482 1.5311 1.4289 1.4125 1 1.3335];
KGM_rho_Maurer = [0.0100 4.2170e-4 5.2723e-05 9.5499e-05 1.4191e-04 100 100 100 1.6141e-04 1.4191e-04];
KGM_m_Maurer = [0 0 0 0 0 0 0 0 0 0];

% SOE
SOE_n_Maurer = [1 1 1 1 1 1 1 1 1 1];
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimate calib

SDR_b = @(K,m,n,phi,T2ML) K ./ (phi.^m.*T2ML.^n);
Seevers_b = @(K,m,n,phi,T2ML,T2B) K ./ (phi.^m.*((T2ML.^(-1) - T2B.^(-1)).^(-1)).^n);
SOE_b = @(K,n,SOE) K ./ (SOE.^n);

for kk = 1:length(sites)
    [T2dist, T2logbins, nmrName] = loadRawNMRdata(sites{kk});
   
    [d, K, T2ML, phi, z, SumEch, logK, logT2ML, logPhi, SumEch_3s, SumEch_twm, ...
    SumEch_twm_3s] = loadnmrdata2(nmrName);
    
    Kall{kk} = K;
    Phiall{kk} = phi;
    T2MLall{kk} = T2ML;
    
    SDR_bFastm0{kk} = SDR_b(K,0,2,phi,T2ML);
    SDR_bSlowm0{kk} = SDR_b(K,0,1,phi,T2ML);
    SDR_bFastm4{kk} = SDR_b(K,4,2,phi,T2ML);
    SDR_bSlowm4{kk} = SDR_b(K,4,1,phi,T2ML);
    SDR_bFastm2{kk} = SDR_b(K,2,2,phi,T2ML);
    SDR_bSlowm2{kk} = SDR_b(K,2,1,phi,T2ML);

    
    %Seevers_bEst{kk} = Seevers_b(K,Seevers_m(kk),Seevers_n(kk),phi,T2ML,T2Bavg);
    SOE_bFast{kk} = SOE_b(K,2,SumEch);
    SOE_bSlow{kk} = SOE_b(K,1,SumEch);

    SumEchAll{kk} = SumEch;
    
end

for kk = 1:length(sites_Maurer)
    [T2dist, T2logbins, nmrName] = loadRawNMRdata(sites_Maurer{kk});
   
    [d, K_Maurer, T2ML_Maurer, phi_Maurer, z, SumEch_Maurer, logK, logT2ML, logPhi, SumEch_3s, SumEch_twm, ...
    SumEch_twm_3s] = loadnmrdata2(nmrName);
    
    Kall_Maurer{kk} = K_Maurer;
    Phiall_Maurer{kk} = phi_Maurer;
    T2MLall_Maurer{kk} = T2ML_Maurer;
    %%%%%%%%%%%%%%%%%%%%
    
    SDR_bFastm0_Maurer{kk} = SDR_b(K_Maurer,0,2,phi_Maurer,T2ML_Maurer);
    SDR_bSlowm0_Maurer{kk} = SDR_b(K_Maurer,0,1,phi_Maurer,T2ML_Maurer);
    SDR_bFastm4_Maurer{kk} = SDR_b(K_Maurer,4,2,phi_Maurer,T2ML_Maurer);
    SDR_bSlowm4_Maurer{kk} = SDR_b(K_Maurer,4,1,phi_Maurer,T2ML_Maurer);
    SDR_bFastm2_Maurer{kk} = SDR_b(K_Maurer,2,2,phi_Maurer,T2ML_Maurer);
    SDR_bSlowm2_Maurer{kk} = SDR_b(K_Maurer,2,1,phi_Maurer,T2ML_Maurer);
    
    %Seevers_bEst{kk} = Seevers_b(K,Seevers_m(kk),Seevers_n(kk),phi,T2ML,T2Bavg);
    SOE_bFast_Maurer{kk} = SOE_b(K_Maurer,2,SumEch_Maurer);
    SOE_bSlow_Maurer{kk} = SOE_b(K_Maurer,1,SumEch_Maurer);

    SumEchAll_Maurer{kk} = SumEch_Maurer;
        
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SDR Slow/Fast Comparision
figure(1)
box on
grid on
hold on

scatter(vertcat(Kall{:}),vertcat(SDR_bFastm0{:}),40,'Filled')
scatter(vertcat(Kall{:}),vertcat(SDR_bSlowm0{:}),40,'Filled')
scatter(vertcat(Kall{:}),vertcat(Kall{:}),20,'Filled')

set(gca,'XScale','log')
set(gca,'YScale','log')
set(gca,'FontSize',16)

xlim([10^-7,10^-2])
ylim([10^-6,10^1])

xlabel('Hydraulic Conductivity (m/s)')
ylabel('Calibration parameter')
title('SDR Slow vs Fast Diffusion Calibration Comparision')
legend({'Wisc SDR n=2, m=0','Wisc SDR n=1, m=0','1:1'},'Location','northwest')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SOE "Slow/Fast" Comparision
figure(2)
box on
grid on
hold on

scatter(vertcat(Kall{:}),vertcat(SOE_bFast{:}),40,'Filled')
scatter(vertcat(Kall{:}),vertcat(SOE_bSlow{:}),40,'Filled')
scatter(vertcat(Kall{:}),vertcat(Kall{:}),20,'Filled')


set(gca,'XScale','log')
set(gca,'YScale','log')
set(gca,'FontSize',16)

xlim([10^-7,10^-2])
ylim([10^-6,10^1])

xlabel('Hydraulic Conductivity (m/s)')
ylabel('Calibration parameter')
legend({'SOE n=2','SOE n=1'},'Location','northwest')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SDR Fast m Comparision
figure(3)
box on
grid on
hold on

scatter(vertcat(Kall{:}),vertcat(SDR_bFastm4{:}),40,'Filled')
scatter(vertcat(Kall{:}),vertcat(SDR_bFastm2{:}),40,'Filled')
scatter(vertcat(Kall{:}),vertcat(SDR_bFastm0{:}),40,'Filled')
scatter(vertcat(Kall{:}),vertcat(Kall{:}),20,'Filled')


set(gca,'XScale','log')
set(gca,'YScale','log')
set(gca,'FontSize',16)

xlim([10^-7,10^-2])
ylim([10^-6,10^1])

xlabel('Hydraulic Conductivity (m/s)')
ylabel('Calibration parameter')
title('Comparing different SDR m, fast diffusion')
legend({'Wisc SDR n=2, m=4','Wisc SDR n=2, m=2','Wisc SDR n=2, m=0','1:1'},'Location','northwest')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SDR Slow m Comparision
figure(4)
box on
grid on
hold on

scatter(vertcat(Kall{:}),vertcat(SDR_bSlowm4{:}),40,'Filled')
scatter(vertcat(Kall{:}),vertcat(SDR_bSlowm2{:}),40,'Filled')
scatter(vertcat(Kall{:}),vertcat(SDR_bSlowm0{:}),40,'Filled')
scatter(vertcat(Kall{:}),vertcat(Kall{:}),20,'Filled')


set(gca,'XScale','log')
set(gca,'YScale','log')
set(gca,'FontSize',16)

xlim([10^-7,10^-2])
ylim([10^-6,10^1])

xlabel('Hydraulic Conductivity (m/s)')
ylabel('Calibration parameter')
legend({'SDR n=1, m=4','SDR n=1, m=2','SDR n=1, m=0'},'Location','northwest')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SDR sites comparision
%%
figure(5)
box on
grid on
hold on

scatter(vertcat(Kall{:}),vertcat(SDR_bFastm0{:}),40,'Filled')
scatter(vertcat(Kall_Maurer{:}),vertcat(SDR_bFastm0_Maurer{:}),40,'Filled')
scatter(vertcat(Kall{:}),vertcat(Kall{:}),20,'Filled')


set(gca,'XScale','log')
set(gca,'YScale','log')
set(gca,'FontSize',16)

xlim([10^-7,10^-2])
ylim([10^-6,10^1])

xlabel('Hydraulic Conductivity (m/s)')
ylabel('Calibration parameter')
title('SDR Calib Comparision')
legend({'Wisc SDR n=2, m=0','K+W SDR n=2, m=0','1:1'},'Location','northwest')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% SOE sites comparision
figure(6)
box on
grid on
hold on

scatter(vertcat(Kall{:}),vertcat(SOE_bSlow{:}),40,'Filled')
scatter(vertcat(Kall_Maurer{:}),vertcat(SOE_bSlow_Maurer{:}),40,'Filled')
scatter(vertcat(Kall{:}),vertcat(Kall{:}),20,'Filled')


set(gca,'XScale','log')
set(gca,'YScale','log')
set(gca,'FontSize',16)

xlim([10^-7,10^-2])
ylim([10^-6,10^1])

xlabel('Hydraulic Conductivity (m/s)')
ylabel('Calibration parameter')
title('SOE Calib Comparision')
legend({'SOE n=1 Wisc','SOE n=1 Maurer','1:1'},'Location','northwest')


%%
% SOE vs K
figure(7)
box on
grid on
hold on

scatter(vertcat(Kall{:}),vertcat(SumEchAll{:}),40,'Filled')
scatter(vertcat(Kall_Maurer{:}),vertcat(SumEchAll_Maurer{:}),40,'Filled')
%scatter(vertcat(Kall{:}),vertcat(Kall{:}),20,'Filled')

set(gca,'XScale','log')
set(gca,'YScale','log')
set(gca,'FontSize',16)

%xlim([10^-7,10^-2])
%ylim([10^-3,10^0])

xlabel('Hydraulic Conductivity (m/s)')
ylabel('SOE')

legend({'SOE n=1 Wisc','SOE n=1 Maurer'},'Location','northwest')
%%%%
%%
% T2ML vs K



T2B = 2.2;
T2ML_Seevers_Wisc = (((1./vertcat(T2MLall{:})) - 1/T2B).^(-1)).^2;
T2ML_Seevers_Maurer = (((1./vertcat(T2MLall_Maurer{:})) + 1/T2B).^(-1)).^2;

Kall_vec = vertcat(Kall{:});
Kall_Maurer_vec = vertcat(Kall_Maurer{:});
T2MLall_vec = vertcat(T2MLall{:});
T2MLall_Maurer_vec = vertcat(T2MLall_Maurer{:});

logKall_vec = log10(Kall_vec);
logKall_Maurer_vec = log10(Kall_Maurer_vec);
logT2MLall_vec = log10(T2MLall_vec.^2);
logT2MLall_Maurer_vec = log10(T2MLall_Maurer_vec.^2);
logT2ML_Seevers_Maurer_vec = log10(T2ML_Seevers_Maurer);

idealT2ML_Wisc = sqrt(Kall_vec./0.0043);
idealT2ML_Maurer = sqrt(Kall_Maurer_vec./0.0278);

% Try to come up with better fit
nPoly = 1;
[pVal, S] = polyfit(logKall_vec,logT2MLall_vec,nPoly);
[pVal_Maurer, S_Maurer] = polyfit(logKall_Maurer_vec,logT2MLall_Maurer_vec,nPoly);

f = fit(logKall_Maurer_vec, logT2MLall_Maurer_vec,'exp1');

% figure(10)
% hold on
% box on 
% grid on
% plot(f,logKall_Maurer_vec,logT2MLall_Maurer_vec)
% 

bestPolyFitT2ML = zeros(length(vertcat(Kall{:})),1);
bestPolyFitT2ML_Maurer = zeros(length(vertcat(Kall_Maurer{:})),1);


for kk = 0:length(pVal)-1
   powerIndex = (length(pVal)-1)-kk;
   polyFitTerm = logKall_vec.^(powerIndex).*pVal(kk+1);
   polyFitTermMaurer = logKall_Maurer_vec.^(powerIndex).*pVal_Maurer(kk+1);
   
   bestPolyFitT2ML = bestPolyFitT2ML + polyFitTerm;
   bestPolyFitT2ML_Maurer = bestPolyFitT2ML_Maurer + polyFitTermMaurer;
      
    %bestPolyFitT2ML = vertcat(Kall{:}).^2*pVal(1)+vertcat(Kall{:})*pVal(2)+pVal(3);
    %bestPolyFitT2ML_Maurer = vertcat(Kall_Maurer{:}).^2*pVal_Maurer(1)+vertcat(Kall_Maurer{:})...
    %    *pVal_Maurer(2)+pVal_Maurer(3);

end

%T2ML_fit_lin = (Kall_Maurer_vec./(10^pVal_Maurer(2))).^(1/(pVal_Maurer(1)));
T2ML_fit_lin = ((Kall_Maurer_vec).^(pVal_Maurer(1))).*10.^(pVal_Maurer(2));




figure(8)
box on
grid on
hold on

%scatter(vertcat(Kall{:}),vertcat(T2MLall{:}).^2,40,'Filled')
scatter(vertcat(Kall_Maurer{:}),vertcat(T2MLall_Maurer{:}).^2,40,'Filled')
%scatter(vertcat(Kall{:}),vertcat(Kall{:}),20,'Filled')

%scatter(vertcat(Kall{:}),(T2ML_Seevers_Wisc(:)),40,'Filled')
%scatter(vertcat(Kall_Maurer{:}),(T2ML_Seevers_Maurer(:)),40,'Filled')

test1 = (idealT2ML_Wisc);
test2 = (idealT2ML_Maurer);

%scatter(vertcat(Kall{:}),test1.^2,20,'Filled')
scatter(vertcat(Kall_Maurer{:}),test2.^2,20,'Filled')

%scatter(vertcat(Kall{:}),10.^bestPolyFitT2ML,20,'r','Filled')
scatter(vertcat(Kall_Maurer{:}),10.^bestPolyFitT2ML_Maurer,20,'r','Filled')
scatter(vertcat(Kall_Maurer{:}),T2ML_fit_lin,20,'b','Filled')


set(gca,'XScale','log')
set(gca,'YScale','log')
set(gca,'FontSize',16)

%xlim([10^-7,10^-2])
%ylim([10^-3,10^0])

xlabel('Hydraulic Conductivity (m/s)')
ylabel('T2ML')

legend({'T2ML^2 Wisc','T2ML^2 Maurer'},'Location','northwest')


%%
% Plot K as a function of T2ML

Kmodel_logFit = 3.58578*(T2MLall_Maurer_vec).^(2/0.48);
Kmodel_SDR_fast = 0.0278.*(T2MLall_Maurer_vec).^(2);
Kmodel_SDR_mod = 10000.*(T2MLall_Maurer_vec).^(6);

Kmodel_SDR_slow = 0.0025.*(T2MLall_Maurer_vec).^(1);

% Compute K diff Factor from new model
kDiffFactorNew = estimateKdiffFactor(Kall_Maurer_vec,Kmodel_logFit,1);
kDiffFactorReg = estimateKdiffFactor(Kall_Maurer_vec,Kmodel_SDR_fast,1);

figure(10)
box on 
grid on
hold on

scatter(T2MLall_Maurer_vec, Kall_Maurer_vec,20,'Filled')

%scatter(T2MLall_Maurer_vec, Kmodel_logFit,20,'Filled')
scatter(T2MLall_Maurer_vec, Kmodel_SDR_fast, 20,'Filled')
scatter(T2MLall_Maurer_vec, Kmodel_SDR_slow, 20,'Filled')
scatter(T2MLall_Maurer_vec, Kmodel_SDR_mod, 20,'Filled')


set(gca,'XScale','log')
set(gca,'YScale','log')
set(gca,'FontSize',16)

xlabel('T2ML')
ylabel('K')


%%

figure(9)
box on
grid on
hold on

%scatter(vertcat(Kall{:}),vertcat(T2MLall{:}).^2,40,'Filled')
scatter(vertcat(Kall_Maurer{:}),vertcat(T2MLall_Maurer{:}).^2,40,'Filled')
%scatter(vertcat(Kall{:}),vertcat(Kall{:}),20,'Filled')

%scatter(vertcat(Kall{:}),(T2ML_Seevers_Wisc(:)),40,'Filled')
scatter(vertcat(Kall_Maurer{:}),(T2ML_Seevers_Maurer(:)),40,'Filled')


%scatter(vertcat(Kall{:}),idealT2ML_Wisc.^2,20,'Filled')
scatter(vertcat(Kall_Maurer{:}),10.^(idealT2ML_Maurer).^2,20,'Filled')




%xlim([10^-7,10^-2])
%ylim([10^-3,10^0])

xlabel('Hydraulic Conductivity (m/s)')
ylabel('T2ML')

legend({'T2ML^2 Wisc','T2ML^2 Maurer'},'Location','northwest')


