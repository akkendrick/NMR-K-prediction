load('wisc_bboot_n2_m0.mat')

figure(1)
nBins = 25;

for kk = 1:length(siteList)
    
    subplot(2,2,kk)
    histogram(b_boot_all{kk},nBins)
    grid on
    box on
    title(siteList{kk})
    
    xlim([1e-4, 1.2e-2])
    ylim([0 250])
    
end
