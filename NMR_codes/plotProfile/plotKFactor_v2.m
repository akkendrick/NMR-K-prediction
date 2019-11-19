function plotKFactor_v2(K,k_est,k_avg,k_names)
    
    kDiffFactor = estimateKdiffFactor(K,k_est,0);
    kDiffFactorAvg = estimateKdiffFactor(K,k_avg,0);

   % kUnitFactor = estimateKdiffFactor(K,K);
    
    figure
    hold on
    
    grid on
    box on
    
    
    [nrow, ncol] = size(k_avg);
    for a = 1:ncol
        plot(K,kDiffFactorAvg(:,a),'^','MarkerSize',8)
    end 
    
    legend(k_names,'Location','northwest')

    [nrow, ncol] = size(k_est);
    for a = 1:ncol
        plot(K,kDiffFactor(:,a),'.k','MarkerSize',2,'HandleVisibility','off')
    end 
    
    plot(K,ones(1,length(K)),'k','LineWidth',2,'HandleVisibility','off')
  
    
    xlabel('DPP K (m/s)')
    ylabel('K Difference Factor') 
    
    set(gca,'FontSize',14)
    set(gca,'XScale','log')
    set(gca,'YScale','log')
    

end

