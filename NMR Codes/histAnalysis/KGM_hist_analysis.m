% KGM Hist
load('KGM_btstrp_Wisc_1000.mat')

for kk = 1:length(siteList)
   figure(kk)
   
   currentTau = bootKGM{1,kk}(:,1);
   currentRho = bootKGM{1,kk}(:,2);
   log10Rho = log10(currentRho);
   
   rhoBins =  linspace(-4.5, -1, 20); 
   tauBins = linspace(1,2.5, 40);
   subplot(121)
   histogram(log10Rho, rhoBins,'Normalization','probability')
   title(strcat('log10(\rho) ',siteList{kk}))
   
   subplot(122)
   histogram(currentTau, tauBins,'Normalization','probability')
   title(strcat('Tau ',siteList{kk}))
   
   median(log10Rho(log10Rho < 2))
   median(currentTau)
   
   
    
end