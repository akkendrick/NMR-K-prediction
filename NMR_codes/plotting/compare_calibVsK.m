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
% sites = {'dpnmr_larned_east','dpnmr_larned_lwph','dpnmr_larned_west',...
%    'dpnmrA11','dpnmrA12','dpnmrC1S','dpnmrC1SE','dpnmrC1SW',...
%    'dpnmr_leque_east','dpnmr_leque_west'}
% 
% load SDR_maurer_bestFit_1015_m0_n2.mat
% 
% sites = siteList;
% SDR_b = squeeze(totalbMatrix(1,1,:))';
% SDR_n = n;
% SDR_n = repmat(n,length(siteList),1);
% SDR_m = squeeze(totalmMatrix(1,1,:))';
% 
% load Seevers_maurer_bestFit_1015_m0_n2.mat
% sites = siteList;
% Seevers_b = squeeze(totalbMatrix(1,1,:))';
% Seevers_n = n;
% Seevers_n = repmat(n,length(siteList),1);
% Seevers_m = squeeze(totalmMatrix(1,1,:))';
% 
% % KGM (has been checked, looks good)
% KGM_tau = [1.6788 1.9498 1 1 1.1482 1.5311 1.4289 1.4125 1 1.3335];
% KGM_rho = [0.0100 4.2170e-4 5.2723e-05 9.5499e-05 1.4191e-04 100 100 100 1.6141e-04 1.4191e-04];
% KGM_m = [0 0 0 0 0 0 0 0 0 0];
% 
% % SOE
% SOE_n = [1 1 1 1 1 1 1 1 1 1];
% SOE_b = [0.0052 0.0025 0.0010 0.0033 0.0032 0.0080 0.0077 0.0074 0.0033 0.0013];

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
    
    SDR_bEst{kk} = SDR_b(K,SDR_m(kk),SDR_n(kk),phi,T2ML);
    Seevers_bEst{kk} = Seevers_b(K,Seevers_m(kk),Seevers_n(kk),phi,T2ML,T2Bavg);
    SOE_bEst{kk} = SOE_b(K,SOE_n(kk),SumEch);
    SumEchAll{kk} = SumEch;
end

numModels = 3;
% for kk = 1:length(sites)
%    
%       
%    figure(1)
%    subplot(length(sites),1,kk)
%    
%    box on
%    grid on
%    hold on
%    
%    scatter(Kall{kk},SDR_bEst{kk},60,'Filled')
%    scatter(Kall{kk},Seevers_bEst{kk},60,'Filled')
%    scatter(Kall{kk},SOE_bEst{kk},60,'Filled')
% 
%    set(gca,'XScale','log')
%    set(gca,'YScale','log')
%    set(gca,'FontSize',16)
% 
%    ylim([5*10^-4,10^-1])
%    xlim([5*10^-6,10^-3])
% 
%    title(sites{kk})
%    xlabel('Hydraulic Conductivity (m/s)')
%    ylabel('Calibration parameter')
%    
%     legend({'SDR','Seevers','SOE'},'Location','northwest')
% end

figure(2)
%subplot(2,1,1)

box on
grid on
hold on

scatter(vertcat(Kall{:}), vertcat(SDR_bEst{:}),60,'Filled')
scatter(vertcat(Kall{:}), vertcat(Seevers_bEst{:}),60,'Filled')
scatter(vertcat(Kall{:}), vertcat(SOE_bEst{:}),60,'Filled')

set(gca,'XScale','log')
set(gca,'YScale','log')
set(gca,'FontSize',16)

%ylim([5*10^-4,10^-1])
%xlim([5*10^-6,10^-3])

xlim([10^-7,10^-2])
ylim([10^-6,10^1])

%title('Kansas + Washington')
title('Wisconsin')

xlabel('Hydraulic Conductivity (m/s)')
ylabel('Calibration parameter')

legend({'SDR','Seevers','SOE'},'Location','northwest')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Try fitting the data
SDR_bAll = vertcat(SDR_bEst{:});
Seevers_bAll = vertcat(Seevers_bEst{:});
SOE_bAll = vertcat(SOE_bEst{:});

Kall_vec = vertcat(Kall{:});
[Kall_vec, ind] = sort(Kall_vec);
logKall_vec = log10(Kall_vec);

SOE_bAll = SOE_bAll(ind);
Seevers_bAll = Seevers_bAll(ind);
SDR_bAll = SDR_bAll(ind);

logSOE_bAll = log10(SOE_bAll);
logSeevers_bAll = log10(Seevers_bAll);
logSDR_bAll = log10(SDR_bAll);

XSOE = [ones(length(SOE_bAll),1) logKall_vec];
logSOE_fit = XSOE\logSOE_bAll;

XSeevers = [ones(length(Seevers_bAll),1) logKall_vec];
logSeevers_fit = XSeevers\logSeevers_bAll;

XSDR = [ones(length(SDR_bAll),1) logKall_vec];
logSDR_fit = XSDR\logSDR_bAll;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2)

fittedlogSOE = XSOE*logSOE_fit;
fittedSOE = 10.^fittedlogSOE;

fittedlogSeevers = XSeevers*logSeevers_fit;
fittedSeevers = 10.^fittedlogSeevers;

fittedlogSDR = XSDR*logSDR_fit;
fittedSDR = 10.^fittedlogSDR;

fittedSOE_take2 = 10.^(logSOE_fit(1)).*Kall_vec.^(logSOE_fit(2));

plot(Kall_vec, fittedSOE,'--')
plot(Kall_vec, fittedSeevers, '--')
plot(Kall_vec, fittedSDR,'--')

plot(Kall_vec, fittedSOE_take2,'--')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure()
box on
grid on
hold on

scatter(vertcat(SumEchAll{:}).^2, vertcat(SOE_bEst{:}),60,'filled')

set(gca,'XScale', 'log')
set(gca,'YScale','log')

