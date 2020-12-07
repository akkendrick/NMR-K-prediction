% Plot results of prior estimates of bounding values of b,m,n

clear
%close all

%load 'SDR_bestFit_table_m0_n2.mat'
%load 'Seevers_bestFit_table.mat'
%load 'TC_bestFit_240_table.mat'
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
load SDR_bestFit_0529_m1_n2.mat
SDR_b = squeeze(totalbMatrix(1,1,:))';
SDR_m = [1 1 1 1];
SDR_n = [2 2 2 2];

load Seevers_bestFit_0422_m1_n2_T2BMap.mat
Seevers_m = [1 1 1 1];
Seevers_n = [2 2 2 2];
Seevers_b = squeeze(totalbMatrix(1,1,:))';

load optimalCutoffTable_n2_m1_RMSE_2000.mat
TC_m = [1 1 1 1];
TC_n = [2 2 2 2];
idealCutoffs = [102.0640, 352.0714, 189.5622, 173.5163]; %in ms 
idealCutoffs = idealCutoffs*10^-3; %in s
refCutoffs = [33 33 33 33]*10^-3;

load KGM_btsrp_Wisc_5.5DegC_m1.mat
KGM_tau = [median(bootKGM{1}(:,1)) median(bootKGM{2}(:,1)) median(bootKGM{3}(:,1))...
    median(bootKGM{4}(:,1))];
KGM_rho = [median(bootKGM{1}(:,2)) median(bootKGM{2}(:,2)) median(bootKGM{3}(:,2))...
    median(bootKGM{4}(:,2))];
KGM_m = [1 1 1 1];

% SOE
load SOE_n1_529.mat
SOE_n = mediann;
SOE_b = medianb;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% K Models
SDR_K = @(b,m,n,phi,T2ML) b.*(phi.^m).*(T2ML).^n;
SOE_K = @(b,n,SOE) b.*(SOE).^n;
Seevers_K = @(b,m,n,T2ML,T2B,phi) b.*(phi).^m.*((T2ML.^(-1) - T2B.^(-1)).^(-1)).^n;

Temp = 5.5;  % temperature in degress C 
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

figureson = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    
   [K,z,T2dist,T2logbins,TC_ideal_Kest{kk},bestFitMatrix,totalError,...
        indexQuotient] = computeTCperm_nobtstp(sites{kk},TC_n(kk),TC_m(kk),median(totalcMatrix(kk,:)),idealCutoffs(kk),figureson);
    
   [K,z,T2dist,T2logbins, TC_ref_Kest{kk},bestFitMatrix,totalError,...
        indexQuotient] = computeTCperm_nobtstp(sites{kk},TC_n(kk),TC_m(kk),median(totalcMatrix(kk,:)),refCutoffs(kk),figureson);
           
    Ktot{kk} = K;
       
end

K_all = vertcat(Ktot{:});
k_estimates = [vertcat(SDR_Kest{:}) vertcat(SOE_Kest{:}) vertcat(Seevers_Kest{:}) vertcat(KGM_Kest{:})...
    vertcat(TC_ideal_Kest{:}) vertcat(TC_ref_Kest{:})];
k_names = {'SDR','SOE','Seevers','KGM','TC Ideal','TC Ref'};

kDiffFactorTemp = estimateKdiffFactor(K_all, k_estimates, 1);

%Figure out good mean, find min
kDiffFactorMean{:,kk} = mean(kDiffFactorTemp);
kDiffFactor{:,kk} = kDiffFactorTemp;

meanK = mean(K_all);
medianK = median(K_all);

figure('Renderer', 'painters', 'Position', [10 10 600 500])
plotKestKdpp_v2(K_all,k_estimates(:,1),k_estimates(:,1),k_names(1))
plotKestKdpp_v2(K_all,k_estimates(:,2),k_estimates(:,2),k_names(1:2))
names = [k_names(1) k_names(2)];
title('Comparing models of K')

ylim([10^-6,5*10^-3])
xlim([10^-6,5*10^-3])

legend(names,'Location','northwest')
fileString = strcat('SDR_SOE','.png');
print('-dpng','-r300',fileString)

figure('Renderer', 'painters', 'Position', [10 10 600 500])
names = [k_names(1) k_names(3)];
plotKestKdpp_v2(K_all,k_estimates(:,1),k_estimates(:,1),names)
plotKestKdpp_v2(K_all,k_estimates(:,3),k_estimates(:,3),names)

title('Comparing models of K')

ylim([10^-6,5*10^-3])
xlim([10^-6,5*10^-3])
legend(names,'Location','northwest')
fileString = strcat('SDR_Seevers','.png');
print('-dpng','-r300',fileString)

figure('Renderer', 'painters', 'Position', [10 10 600 500])
names = [k_names(1) k_names(4)];
plotKestKdpp_v2(K_all,k_estimates(:,1),k_estimates(:,1),names)
plotKestKdpp_v2(K_all,k_estimates(:,4),k_estimates(:,4),names)

title('Comparing models of K')

ylim([10^-6,5*10^-3])
xlim([10^-6,5*10^-3])

legend(names,'Location','northwest')
fileString = strcat('SDR_KGM','.png');
print('-dpng','-r300',fileString)

figure('Renderer', 'painters', 'Position', [10 10 600 500])
names = [k_names(1) k_names(5)];
plotKestKdpp_v2(K_all,k_estimates(:,1),k_estimates(:,1),names)
plotKestKdpp_v2(K_all,k_estimates(:,5),k_estimates(:,5),names)

title('Comparing models of K')
meanErrorFactorIdeal = mean(estimateKdiffFactor(K_all,k_estimates(:,5),1));
medianErrorFactorIdeal = median(estimateKdiffFactor(K_all,k_estimates(:,5),1))
errorFactor = estimateKdiffFactor(K_all,k_estimates(:,6),1);

ylim([10^-6,5*10^-3])
xlim([10^-6,5*10^-3])

legend(names,'Location','northwest')
fileString = strcat('SDR_TC_ideal','.png');
print('-dpng','-r300',fileString)

errorFactor(errorFactor > 10) = 10;

figure()
histogram(errorFactor)

figure('Renderer', 'painters', 'Position', [10 10 600 500])
names = [k_names(1) k_names(6)];
plotKestKdpp_v2(K_all,k_estimates(:,1),k_estimates(:,1),names)
plotKestKdpp_v2(K_all,k_estimates(:,6),k_estimates(:,6),names)

title('Comparing models of K')
meanErrorFactorRef = mean(estimateKdiffFactor(K_all,k_estimates(:,6),1));
medianErrorFactorRef = median(estimateKdiffFactor(K_all,k_estimates(:,6),1))
errorFactor = estimateKdiffFactor(K_all,k_estimates(:,6),1);

ylim([10^-6,5*10^-3])
xlim([10^-6,5*10^-3])

legend(names,'Location','northwest')
fileString = strcat('SDR_TC_ref','.png');
print('-dpng','-r300',fileString)

errorFactor(errorFactor > 10) = 10;

figure()
histogram(errorFactor)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
