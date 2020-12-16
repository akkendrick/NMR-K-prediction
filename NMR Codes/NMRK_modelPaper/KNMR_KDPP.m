% Make NMR K vs  DPP K plot
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
%TC model

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
%

%%
totalSDRK = vertcat(SDRK{:});
totalSeeversK = vertcat(SeeversK{:});
totalKGMK = vertcat(KGMK{:});
totalSOEK = vertcat(SOEK{:});
totalTCK = vertcat(totalkTC{:});
totalTCK = vertcat(totalTCK{:});

totalDPPK = vertcat(DPP_K{:});

figure(1)

hold on
    
grid on
box on

scatter(totalDPPK, totalSDRK,60,'Filled')
scatter(totalDPPK, totalSeeversK,60,'Filled')
scatter(totalDPPK, totalKGMK,60,'Filled')
scatter(totalDPPK, totalSOEK,60,'Filled')
scatter(totalDPPK, totalTCK,60,'Filled')

legend('SDR','Seevers','KGM','SOE','TC','Location','northwest')

plot(totalDPPK,totalDPPK,'k','LineWidth',2,'HandleVisibility','off')
plot(totalDPPK,totalDPPK*10,'k:','HandleVisibility','off')
plot(totalDPPK,totalDPPK*0.1,'k:','HandleVisibility','off')

xlabel('DPP K (m/s)')
ylabel('Estimated K (m/s)') 

set(gca,'FontSize',16)
set(gca,'XScale','log')
set(gca,'YScale','log')
