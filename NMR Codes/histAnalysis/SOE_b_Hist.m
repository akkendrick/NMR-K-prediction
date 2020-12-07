load('SOE_n1_Wisc.mat')

figure(1)
nBins = 25;

for kk = 1:length(siteList)
    
    subplot(2,2,kk)
    histogram(b_boot_all{kk},nBins)
    grid on
    box on
    title(siteList{kk})
    
    xlim([1e-4, 2.5e-2])
    ylim([0 150])
    
end
