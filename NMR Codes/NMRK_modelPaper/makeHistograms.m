% Make nice histograms of 4 NMR-K models on Wisconsin data
clear
close all

%12/7/20 
%Using m = 1, n = 2 IF THIS IS CHANGED CHECK MAIN MODEL FILE CODE, NEED TO
%UPDATE VALUES IN SCRIPTS
m = 1;
n = 2;
siteList = [{"Site1-WellG5"},{"Site1-WellG6"},{"Site2-WellPN1"},{"Site2-WellPN2"}];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% SDR model

runModelPairs

SDRm = [m m m m]';
SDRn = [n n n n]';
SDRb = squeeze(totalbMatrix(1,1,:));
SDRK = totalKMatrix;
DPP_K = DPP_KMatrix;

for k = 1:length(siteList)
  
   [SDR_errorSign,SDR_errorFactor] = estimateKdiffFactor_withSign(DPP_K{:,:,k},SDRK{:,:,k},1);
   SDR_errorSigns{:,k} = SDR_errorSign;
   SDR_errorFactors{:,k} = SDR_errorFactor;
 
end
%%
% SOE Model

computeSOE

SOEn = mediann;
SOEb = medianb;
SOEK = k_estimates;
%%
for k = 1:length(siteList)
   k
   
   [SOE_errorSign,SOE_errorFactor] = estimateKdiffFactor_withSign(DPP_K{:,:,k},SOEK{:,k},1);
   SOE_errorSigns{:,k} = SOE_errorSign;
   SOE_errorFactors{:,k} = SOE_errorFactor;
 
end

%%
% KGM Model

calcKGM

KGMK = KGM_K;
KGM_errorSigns = errorSign;
KGM_errorFactors = errorFactor;

%%
% Seevers Model 

Seeversm = [m m m m]';
Seeversn = [n n n n]';

for k = 1:length(siteList)
    
    site = siteList{k}
    [K,z,T2dist,T2logbins,SeeversK{k},k_mcmc,k_direct,bestFitMatrix,b_boot,totalErrorEstimate] = computeSeevers(site,Seeversn(k),Seeversm(k),figureson,wDirect);
   
   [Seevers_errorSign,Seevers_errorFactor] = estimateKdiffFactor_withSign(DPP_K{:,:,k},SeeversK{k},1);
   Seevers_errorSigns{:,k} = Seevers_errorSign;
   Seevers_errorFactors{:,k} = Seevers_errorFactor;
end

%%
% %TC model 42.1
% 
% TCm = [m m m m];
% TCn = [n n n n];
% 
% TC_gridSearch_wisc;
% 
% cutoff = [cutoff cutoff cutoff cutoff];
% cTC = totalcMatrix;
% 
% for k = 1:length(siteList)
% 
%     [TC_errorSign,TC_errorFactor] = estimateKdiffFactor_withSign(DPP_K{k},totalkTC{k}{:},1);
%     TC_errorSigns{:,k} = TC_errorSign;
%     TC_errorFactors{:,k} = TC_errorFactor;
%     
% end
% 
% %TC model 33
% 
% TCm = [m m m m];
% TCn = [n n n n];
% 
% TC_gridSearch_wisc_33;
% 
% cutoffRef = [cutoff cutoff cutoff cutoff];
% cTCRef = totalcMatrix;
% 
% for k = 1:length(siteList)
% 
%     [TC_errorSign,TC_errorFactor] = estimateKdiffFactor_withSign(DPP_K{k},totalkTC{k}{:},1);
%     TC_errorSigns_ref{:,k} = TC_errorSign;
%     TC_errorFactors_ref{:,k} = TC_errorFactor;
%     
% end
% %

%%
% Check if higher errors are present
% plottedTCerrorFactor = vertcat(TC_errorFactors{:});
% plottedTCerrorFactor_ref = vertcat(TC_errorFactors_ref{:});
plottedSDRerrorFactor = vertcat(SDR_errorFactors{:});
plottedSeeverserrorFactor = vertcat(Seevers_errorFactors{:});
plottedKGMerrorFactor = vertcat(KGM_errorFactors{:});
plottedSOEerrorFactor = vertcat(SOE_errorFactors{:});

% maxTCerrorFactor = max(plottedTCerrorFactor)
% maxTCerrorFactor_ref = max(plottedTCerrorFactor_ref)
medianSDRerrorFactor = median(plottedSDRerrorFactor)
medianSOEerrorFactor = median(plottedSOEerrorFactor)
medianKGMerrorFactor = median(plottedKGMerrorFactor)
medianSeeverserrorFactor = median(plottedSeeverserrorFactor)

%%
% Fix higher errors to make sure they show up on the plot, set +/-;

% plottedTCerrorFactor(plottedTCerrorFactor >= 100) = 200;
% %plottedTCerrorFactor = plottedTCerrorFactor.*vertcat(TC_errorSigns{:});
% 
% plottedTCerrorFactor_ref(plottedTCerrorFactor_ref >= 100) = 200;

plottedSDRerrorFactor(plottedSDRerrorFactor >= 100) = 60;
%plottedSDRerrorFactor = plottedSDRerrorFactor.*vertcat(SDR_errorSigns{:});

plottedSeeverserrorFactor(plottedSeeverserrorFactor >= 100) = 60;
%plottedSeeverserrorFactor = plottedSeeverserrorFactor.*vertcat(Seevers_errorSigns{:});

plottedKGMerrorFactor(plottedKGMerrorFactor >= 100) = 60;
%plottedKGMerrorFactor = plottedKGMerrorFactor.*vertcat(KGM_errorSigns{:});

plottedSOEerrorFactor(plottedSOEerrorFactor >= 100) = 60;
%plottedSOEerrorFactor = plottedSOEerrorFactor.*vertcat(SOE_errorSigns{:});

figure(1)
subplot(2,2,1)
title('SDR Model')

edgesLog = [-2 -1.5 -1 -0.8 -0.6 -0.4 -0.2 0 0.2 0.4 0.6 0.8 1 1.5 2];
edges = 10.^edgesLog;

hold on
grid on
box on
grid minor

histDataSDR = log10(plottedSDRerrorFactor);
histDataSDR = histDataSDR.*vertcat(SDR_errorSigns{:});
histogram(histDataSDR,edgesLog)
ylim([0 50])
%ylim([0 30])

text(-1.9,40,strcat('Median K_{diff} = ',num2str(round(medianSDRerrorFactor))),'FontSize',14)
xticks([-2,-1,0,1,2])%in log space
xticklabels({'-100','-10','0','10','100'})

xlabel('K Difference Factor')
ylabel('Counts')
set(gca,'FontSize',14)

subplot(2,2,2)
title('SOE Model')

edgesLog = [-2 -1.5 -1 -0.8 -0.6 -0.4 -0.2 0 0.2 0.4 0.6 0.8 1 1.5 2];
edges = 10.^edgesLog;

hold on
grid on
box on
grid minor
text(-1.9,40,strcat('Median K_{diff} = ',num2str(round(medianSOEerrorFactor))),'FontSize',14)
histDataSOE = log10(plottedSOEerrorFactor);
histDataSOE = histDataSOE.*vertcat(SOE_errorSigns{:});
histogram(histDataSOE,edgesLog)
ylim([0 50])
%ylim([0 30])
xticks([-2,-1,0,1,2])%in log space

xticklabels({'-100','-10','0','10','100'})
xlabel('K Difference Factor')
ylabel('Counts')
set(gca,'FontSize',14)

subplot(2,2,3)
title('Seevers Model')

edgesLog = [-2 -1.5 -1 -0.8 -0.6 -0.4 -0.2 0 0.2 0.4 0.6 0.8 1 1.5 2];
edges = 10.^edgesLog;

hold on
grid on
box on
grid minor

histDataSeevers = log10(plottedSeeverserrorFactor);
histDataSeevers = histDataSeevers.*vertcat(Seevers_errorSigns{:});
histogram(histDataSeevers,edgesLog)
ylim([0 50])
%ylim([0 30])
xticks([-2,-1,0,1,2])%in log space
xticklabels({'-100','-10','0','10','100'})
xlabel('K Difference Factor')
ylabel('Counts')
set(gca,'FontSize',14)
text(-1.9,40,strcat('Median K_{diff} = ',num2str(round(medianSeeverserrorFactor))),'FontSize',14)

subplot(2,2,4)
title('KGM Model')

edgesLog = [-2 -1.5 -1 -0.8 -0.6 -0.4 -0.2 0 0.2 0.4 0.6 0.8 1 1.5 2];
edges = 10.^edgesLog;

hold on
grid on
box on
grid minor

histDataKGM = log10(plottedKGMerrorFactor);
histDataKGM = histDataKGM.*vertcat(KGM_errorSigns{:});
histogram(histDataKGM,edgesLog)
ylim([0 50])
%ylim([0 30])
text(-1.9,40,strcat('Median K_{diff} = ',num2str(round(medianKGMerrorFactor))),'FontSize',14)

xticks([-2,-1,0,1,2])%in log space
xticklabels({'-100','-10','0','10','100'})
xlabel('K Difference Factor')
ylabel('Counts')
set(gca,'FontSize',14)

% subplot(2,3,5)
% title('TC Model 42.1 ms cutoff')
% 
% edgesLog = [-2 -1.5 -1 -0.8 -0.6 -0.4 -0.2 0 0.2 0.4 0.6 0.8 1 1.5 2];
% edges = 10.^edgesLog;
% 
% hold on
% grid on
% box on
% grid minor
% 
% histDataTC = log10(plottedTCerrorFactor);
% histDataTC = histDataTC.*vertcat(TC_errorSigns{:});
% histogram(histDataTC,edgesLog)
% ylim([0 50])
% %ylim([0 30])
% text(-1.9,40,strcat('Max K_{diff} = ',num2str(round(maxTCerrorFactor))),'FontSize',14)
% 
% xticks([-2,-1,0,1,2])%in log space
% xticklabels({'-100','-10','0','10','100'})
% xlabel('K Difference Factor')
% ylabel('Counts')
% set(gca,'FontSize',14)
% 
% subplot(2,3,6)
% title('TC Model 33 ms cutoff')
% 
% edgesLog = [-2 -1.5 -1 -0.8 -0.6 -0.4 -0.2 0 0.2 0.4 0.6 0.8 1 1.5 2];
% edges = 10.^edgesLog;
% 
% hold on
% grid on
% box on
% grid minor
% 
% histDataTC = log10(plottedTCerrorFactor_ref);
% histDataTC = histDataTC.*vertcat(TC_errorSigns_ref{:});
% histogram(histDataTC,edgesLog)
% ylim([0 50])
% %ylim([0 30])
% text(-1.9,40,strcat('Max K_{diff} = ',num2str(round(maxTCerrorFactor_ref))),'FontSize',14)
% 
% xticks([-2,-1,0,1,2])%in log space
% xticklabels({'-100','-10','0','10','100'})
% xlabel('K Difference Factor')
% ylabel('Counts')
% set(gca,'FontSize',14)
