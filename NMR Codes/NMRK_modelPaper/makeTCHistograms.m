% Make nice histograms of TC model on Wisconsin data with different cutoff
% times
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
%TC model 42.1

TCm = [m m m m];
TCn = [n n n n];

TC_gridSearch_wisc;

cutoff = [cutoff cutoff cutoff cutoff];
cTC = totalcMatrix;

for k = 1:length(siteList)

    [TC_errorSign,TC_errorFactor] = estimateKdiffFactor_withSign(DPP_K{k},totalkTC{k}{:},1);
    TC_errorSigns{:,k} = TC_errorSign;
    TC_errorFactors{:,k} = TC_errorFactor;
    
end
%%
%TC model 33

TCm = [m m m m];
TCn = [n n n n];

TC_gridSearch_wisc_33;

cutoffRef = [cutoff cutoff cutoff cutoff];
cTCRef = totalcMatrix;

for k = 1:length(siteList)

    [TC_errorSign,TC_errorFactor] = estimateKdiffFactor_withSign(DPP_K{k},totalkTC{k}{:},1);
    TC_errorSigns_ref{:,k} = TC_errorSign;
    TC_errorFactors_ref{:,k} = TC_errorFactor;
    
end
%%
%TC model opt

TCm = [m m m m];
TCn = [n n n n];
%% THIS IS SLOW, DON'T RUN UNLESS NECESSARY
clear totalkTC
TC_gridSearch_wisc_opt;

% %Fast version
% load TC_opt_out
%%
cutoffOpt = [cutoff cutoff cutoff cutoff];
cTCOpt = totalcMatrix;

for k = 1:length(siteList)

    [TC_errorSign,TC_errorFactor] = estimateKdiffFactor_withSign(DPP_K{k},bestKEstimate{k},1);
    TC_errorSigns_opt{:,k} = TC_errorSign;
    TC_errorFactors_opt{:,k} = TC_errorFactor;
    
    mean(TC_errorFactors_opt{:,k})
    median(TC_errorFactors_opt{:,k})
    
end

%

%%
% Check if higher errors are present
plottedTCerrorFactor_opt = vertcat(TC_errorFactors_opt{:});
plottedTCerrorFactor = vertcat(TC_errorFactors{:});
plottedTCerrorFactor_ref = vertcat(TC_errorFactors_ref{:});
plottedSDRerrorFactor = vertcat(SDR_errorFactors{:});

medianTCerrorFactor = median(plottedTCerrorFactor);
medianTCerrorFactor_ref = median(plottedTCerrorFactor_ref);
medianTCerrorFactor_opt = median(plottedTCerrorFactor_opt);
medianSDRerrorFactor = median(plottedSDRerrorFactor);


%%
% Fix higher errors to make sure they show up on the plot, set +/-;

plottedTCerrorFactor(plottedTCerrorFactor >= 100) = 60;
%plottedTCerrorFactor = plottedTCerrorFactor.*vertcat(TC_errorSigns{:});

plottedTCerrorFactor_ref(plottedTCerrorFactor_ref >= 100) = 60;

plottedTCerrorFactor_opt(plottedTCerrorFactor_opt >= 100) = 60;

plottedSDRerrorFactor(plottedSDRerrorFactor >= 100) = 60;
%plottedSDRerrorFactor = plottedSDRerrorFactor.*vertcat(SDR_errorSigns{:});


figure()
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

text(-1.2,40,strcat('Median K_{diff} = ',num2str(round(medianSDRerrorFactor))),'FontSize',14)
xticks([-2,-1,0,1,2])%in log space
xticklabels({'-100','-10','0','10','100'})

xlabel('K Difference Factor')
ylabel('Counts')
set(gca,'FontSize',14)

subplot(2,2,4)
title('TC Model 42.1 ms cutoff')

hold on
grid on
box on
grid minor

histDataTC = log10(plottedTCerrorFactor);
histDataTC = histDataTC.*vertcat(TC_errorSigns{:});
histogram(histDataTC,edgesLog)
ylim([0 50])
%ylim([0 30])
text(-1.2,40,strcat('Median K_{diff} = ',num2str(round(medianTCerrorFactor))),'FontSize',14)

xticks([-2,-1,0,1,2])%in log space
xticklabels({'-100','-10','0','10','100'})
xlabel('K Difference Factor')
ylabel('Counts')
set(gca,'FontSize',14)

subplot(2,2,3)
title('TC Model 33 ms cutoff')

edgesLog = [-2 -1.5 -1 -0.8 -0.6 -0.4 -0.2 0 0.2 0.4 0.6 0.8 1 1.5 2];
edges = 10.^edgesLog;

hold on
grid on
box on
grid minor

histDataTC = log10(plottedTCerrorFactor_ref);
histDataTC = histDataTC.*vertcat(TC_errorSigns_ref{:});
histogram(histDataTC,edgesLog)
ylim([0 50])
%ylim([0 30])
text(-1.2,40,strcat('Median K_{diff} = ',num2str(round(medianTCerrorFactor_ref))),'FontSize',14)

xticks([-2,-1,0,1,2])%in log space
xticklabels({'-100','-10','0','10','100'})
xlabel('K Difference Factor')
ylabel('Counts')
set(gca,'FontSize',14)

subplot(2,2,2)
title('TC Model Optimal Cutoff')

edgesLog = [-2 -1.5 -1 -0.8 -0.6 -0.4 -0.2 0 0.2 0.4 0.6 0.8 1 1.5 2];
edges = 10.^edgesLog;

hold on
grid on
box on
grid minor

histDataTC = log10(plottedTCerrorFactor_opt);
histDataTC = histDataTC.*vertcat(TC_errorSigns_opt{:});
histogram(histDataTC,edgesLog)
ylim([0 50])
%ylim([0 30])
text(-1.2,40,strcat('Median K_{diff} = ',num2str(round(medianTCerrorFactor_opt))),'FontSize',14)

xticks([-2,-1,0,1,2])%in log space
xticklabels({'-100','-10','0','10','100'})
xlabel('K Difference Factor')
ylabel('Counts')
set(gca,'FontSize',14)
    