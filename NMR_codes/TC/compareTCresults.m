% Compare different Timur Coates cutoff estimates 
% Figure out which combination of n, m give a cutoff that minimizes error
% the best


data1 = load('optimalCutoffTable_n1_m0.mat'); 
data2 = load('optimalCutoffTable_n1_m2.mat'); 
data3 = load('optimalCutoffTable_n1_m4.mat'); 



% Plot the data and compare

for jj = 1:length(data1.siteList)
    
    figure(1)
    subplot(2,2,jj)
    
    hold on
    
    plot(data1.cutoff,data1.totalErrorMatrix(jj,:))
    plot(data2.cutoff,data2.totalErrorMatrix(jj,:))
    plot(data3.cutoff,data3.totalErrorMatrix(jj,:))

    title(data1.siteList{jj})
    
    ylim([10^-5,10^-2])
    xlim([0, 1])
    
    xlabel('Cutoff (ms)')
    ylabel('Error')
    
    grid on
    box on
    
    set(gca,'FontSize',14)
    
    
end