function plotKFactor(K,k_estimates,k_names,k_sym)
    
    kDiffFactor = estimateKdiffFactor(K,k_estimates,0);
   % kUnitFactor = estimateKdiffFactor(K,K);
    
    figure
    hold on
    
    grid on
    box on
    
    plot(K,ones(1,length(K)),'LineWidth',2)
    
    [nrow, ncol] = size(k_estimates);
    for a = 1:ncol
        plot(K,kDiffFactor(:,a),k_sym{a},'MarkerSize',8)
    end 
    
    legend(k_names,'Location','northwest')
    
    xlabel('DPP K (m/s)')
    ylabel('K Difference Factor') 
    
    set(gca,'FontSize',14)
    set(gca,'XScale','log')
    set(gca,'YScale','log')
    

end

