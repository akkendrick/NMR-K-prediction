% Plot K diff factor for different models, copmare to histogram of data 

clear
%close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sites = [{'Site1-WellG6'} {'Site1-WellG5'} {'Site2-WellPN1'} {'Site2-WellPN2'}];
% figureNames = {'Adams G6','Adams G5','Plainfield PN1','Plainfield PN2'}
% 
% T2B = 2.2293;
% 
% load SDR_bestFit_1101_m1_n2.mat
% sites = siteList;
% SDR_b = squeeze(totalbMatrix(1,1,:))';
% SDR_n = n;
% SDR_n = repmat(n,length(siteList),1);
% SDR_m = squeeze(totalmMatrix(1,1,:))';
% 
% load Seevers_bestFit_1101_m1_n2_T2BAvg.mat
% sites = siteList;
% Seevers_b = squeeze(totalbMatrix(1,1,:))';
% Seevers_n = n;
% Seevers_n = repmat(n,length(siteList),1);
% Seevers_m = squeeze(totalmMatrix(1,1,:))';
% 
% % KGM (has been checked, looks good)
% KGM_tau = [1 1 1.135 1.5668];
% KGM_rho = [4.3251e-05 6.7421e-05 5.6754e-04 1.00e+02];
% KGM_m = [1 1 1 1];
% 
% % SOE
% SOE_n = [1 1 1 1];
% SOE_b = [0.0045 0.0052 0.0158 0.0092];
% 
% Temp = 7.4;  % temperature in degress C 
% Tb = @(Tt) T2B;      % seconds
% 
% plotIndex = 1;
% plotTitle = 'Wisconsin';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sites = {'dpnmr_larned_east','dpnmr_larned_lwph','dpnmr_larned_west',...
   'dpnmrA11','dpnmrA12','dpnmrC1S','dpnmrC1SE','dpnmrC1SW',...
   'dpnmr_leque_east','dpnmr_leque_west'}

figureNames = {'Larned E','Larned lwph','Larned W','A11','A12','C1S','C1SE','C1SW','Leque E','Leque W'}

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

Temp = 12.5;  % temperature in degress C 
Tb = @(Tt) 3.3 + 0.044*(Tt - 35);       % seconds
plotIndex = 2;
plotTitle = 'Kansas + Washington';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% K estimates
SDR_K = @(b,m,n,phi,T2ML) b.*(phi.^m).*(T2ML).^n;
SOE_K = @(b,n,SOE) b.*(SOE).^n;
Seevers_K = @(b,m,n,T2ML,T2B,phi) b.*(phi).^m.*((T2ML.^(-1) - T2B.^(-1)).^(-1)).^n;

rho = @(Tt) 1000*(1 - ((Tt+288.94)./(508929*(Tt+68.12))).*(Tt-3.98).^2); % kg/m^3
eta = @(Tt) 0.0013 - 1.7e-5*Tt;         % Pa -s
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
    
   [SDR_sign{kk}, SDR_diffFactor{kk}] = estimateKdiffFactor_withSign(K,SDR_Kest{kk},1);
   [Seevers_sign{kk}, Seevers_diffFactor{kk}] = estimateKdiffFactor_withSign(K,Seevers_Kest{kk},1);
   [SOE_sign{kk}, SOE_diffFactor{kk}] = estimateKdiffFactor_withSign(K,SOE_Kest{kk},1);
   [KGM_sign{kk}, KGM_diffFactor{kk}] = estimateKdiffFactor_withSign(K,KGM_Kest{kk},1);
    
    SDR_diffFactor{kk}(SDR_diffFactor{kk} > 100) = 100;
    Seevers_diffFactor{kk}(Seevers_diffFactor{kk} > 100) = 100;
    SOE_diffFactor{kk}(SOE_diffFactor{kk} > 100) = 100;
    KGM_diffFactor{kk}(KGM_diffFactor{kk} > 100) = 100;

    SumEch_sites{kk} = SumEch;
    T2ML_sites{kk} = T2ML;
    phi_sites{kk} = phi;
           
end

figure(1)

subplot(2,1,plotIndex)

hold on
grid on 
box on

for kk = 1:length(sites)
    
    plot(kk, mean(SDR_diffFactor{kk}),'k+','MarkerSize',10,'LineWidth',2)
    plot(kk, mean(Seevers_diffFactor{kk}),'rs','MarkerSize',10,'LineWidth',2)
    plot(kk, mean(SOE_diffFactor{kk}),'bo','MarkerSize',10,'LineWidth',2)
    plot(kk, mean(KGM_diffFactor{kk}),'m^','MarkerSize',10,'LineWidth',2)
  
    kk
    
end

%names = {'Adams G6','Adams G5','Plainfield PN1','Plainfield PN2'};

names = figureNames;

xlim([0,length(sites)+1])
ylim([1,100])
set(gca,'xtick',[1:length(sites)],'xticklabel',names)

ylabel('K Difference Factor')
set(gca,'FontSize',16)

legend({'SDR','Seevers','SOE','KGM'})

%mean(vertcat(KGM_diffFactor{:}))

set(gca,'YScale','log')

% Compute the mean and median difference factor across all sites
mean(vertcat(SDR_diffFactor{:}))
mean(vertcat(Seevers_diffFactor{:}))
mean(vertcat(SOE_diffFactor{:}))
mean(vertcat(KGM_diffFactor{:}))

median(vertcat(SDR_diffFactor{:}))
median(vertcat(Seevers_diffFactor{:}))
median(vertcat(SOE_diffFactor{:}))
median(vertcat(KGM_diffFactor{:}))

%
display('Estimating number above factor of 10')
SDR_diffFactor_all = vertcat(SDR_diffFactor{:});
Seevers_diffFactor_all = vertcat(Seevers_diffFactor{:});
SOE_diffFactor_all = vertcat(SOE_diffFactor{:});
KGM_diffFactor_all = vertcat(KGM_diffFactor{:});

SDR_sign_all = vertcat(SDR_sign{:});
Seevers_sign_all = vertcat(Seevers_sign{:});
SOE_sign_all = vertcat(SOE_sign{:});
KGM_sign_all = vertcat(KGM_sign{:});

SumEch_all = vertcat(SumEch_sites{:});
phi_all = vertcat(phi_sites{:});
T2ML_all = vertcat(T2ML_sites{:});

K_all = vertcat(Ktot{:});

display('SDR')
length(find(SDR_diffFactor_all > 10))

display('Seevers')
length(find(Seevers_diffFactor_all > 10))

display('SOE')
length(find(SOE_diffFactor_all > 10))

display('KGM')
length(find(KGM_diffFactor_all > 10))


%%
% Try plotting the difference factors vs K 

figure()
subplot(2,1,1)
hold on
grid on
box on

% scatter(vertcat(Ktot{:}),vertcat(SDR_diffFactor{:}),40,'Filled')
% scatter(vertcat(Ktot{:}),vertcat(Seevers_diffFactor{:}),40,'Filled')
% scatter(vertcat(Ktot{:}),vertcat(SOE_diffFactor{:}),40,'Filled')
% scatter(vertcat(Ktot{:}),vertcat(KGM_diffFactor{:}),40,'Filled')

scatter(vertcat(Ktot{:}),SDR_diffFactor_all.*SDR_sign_all,40,'Filled')
%scatter(vertcat(Ktot{:}),Seevers_diffFactor_all.*Seevers_sign_all,40,'Filled')
scatter(vertcat(Ktot{:}),SOE_diffFactor_all.*SOE_sign_all,40,'Filled')
%scatter(vertcat(Ktot{:}),KGM_diffFactor_all.*KGM_sign_all,40,'Filled')

ylim([-20, 20])


set(gca,'XScale','log')
%set(gca,'YScale','log')

%title(plotTitle)
%title('Wisconsin')

xlabel('Hydraulic Conductivity (m/s)')
ylabel('K difference factor')

xlim([10^-7,10^-2])
%xlim([2e-06,10^-3])

set(gca,'FontSize',16)
legend({'SDR','Seevers','SOE','KGM'})

subplot(2,1,2)

%edges = logspace(-6, -2.3, 20);
edges = logspace(-7, -2, 20);

goodData = K_all(SDR_diffFactor_all < 10);

histogram(vertcat(Ktot{:}),edges);
hold on
histogram(goodData, edges);

set(gca,'XScale','log')
grid on
box on

%xlim([10^-6,5*10^-3])
xlim([10^-7,10^-2])
%xlim([2e-06,10^-3])

set(gca,'FontSize',16)
ylabel('Counts')
xlabel('Hydraulic Conductivity (m/s)')

%%
% Try discretizing data into bins
allK = vertcat(Ktot{:});
allLogK = log10(allK);
allSeeversDiff = vertcat(Seevers_diffFactor{:});
allSDRDiff = vertcat(SDR_diffFactor{:});
allSOEDiff = vertcat(SOE_diffFactor{:});
allKGMDiff = vertcat(KGM_diffFactor{:});

Nbins = 20;
[discretK, edges] = discretize(allLogK, Nbins);

for kk = 1:Nbins
   Kcount(kk) = length(allLogK(discretK == kk));
   subsetSeevers{kk} = allSeeversDiff(discretK == kk);
   subsetSDR{kk} = allSDRDiff(discretK == kk);
   subsetSOE{kk} = allSOEDiff(discretK == kk);
   subsetKGM{kk} = allKGMDiff(discretK == kk);
   
   binnedMeanSeeversDiff(kk) = mean(subsetSeevers{kk});
   binnedMeanSDRDiff(kk) = mean(subsetSDR{kk});
   binnedMeanSOEDiff(kk) = mean(subsetSOE{kk});
   binnedMeanKGMDiff(kk) = mean(subsetKGM{kk});
   
end

X = 1:Nbins;

figure(3)
subplot(2,1,plotIndex)
hold on
grid on
box on

scatter(X,binnedMeanSDRDiff, 40, 'Filled')
scatter(X,binnedMeanSeeversDiff, 40, 'Filled')
scatter(X,binnedMeanSOEDiff, 40, 'Filled')
scatter(X,binnedMeanKGMDiff, 40, 'Filled')

set(gca,'YScale','log')
set(gca,'FontSize',14)
legend({'SDR','Seevers','SOE','KGM'})
xticks(X)
xticklabels(string(edges(1:21)))

ylabel('Mean K Diff Factor')
xlabel('Log(K) Bound')
title(plotTitle)


T2mlPhi = (T2ML_all.^2) .* phi_all;
figure()
scatter(T2mlPhi, SumEch_all, 40,'Filled')
grid on
box on
%set(gca,'XScale','log')
xlabel('(T2ML^2)(\phi) (s)')
ylabel('Sum of Spin-Echoes')

figure()
scatter(T2mlPhi, K_all, 40,'Filled')
grid on
box on
set(gca,'YScale','log')
set(gca,'XScale','log')

xlabel('(T2ML)(\phi) (s)')
ylabel('K (m/s)')

figure()
scatter(SumEch_all, K_all, 40,'Filled')
grid on
box on
set(gca,'YScale','log')
set(gca,'XScale','log')
xlabel('Sum of Spin-Echoes')
ylabel('K (m/s)')


