function plotKestKdpp(K,k_estimates,k_names,k_sym)
    
    figure
    hold on
    
    grid on
    box on
    
    plot(K,K,'LineWidth',2)
    plot(K,K*10,'k:')
    plot(K,K*0.1,'k:')
    
    [nrow, ncol] = size(k_estimates);
    for a = 1:ncol
        plot(K,k_estimates(:,a),k_sym{a},'MarkerSize',8)
    end 
    
    legend(k_names,'Location','northwest')
    
    xlabel('DPP K (m/s)')
    ylabel('Estimated K (m/s)') 
    
    set(gca,'FontSize',14)
    set(gca,'XScale','log')
    set(gca,'YScale','log')
    
%     ylim([10^-6,10^-3])
%     xlim([10^-6,10^-3])
end

