load('Seevers_bestFit_1201_m1_n2_T2BAvg.mat')

figure(1)
nBins = 25;

for kk = 1:length(siteList)
    
    subplot(2,2,kk)
    histogram(b_boot_all{kk},nBins)
    grid on
    box on
    title(siteList{kk})
    
    xlim([1e-4, 3e-3])
    ylim([0 300])
    
end
