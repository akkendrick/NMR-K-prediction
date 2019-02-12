function plotKestKdpp(K,k_estimates,k_names,k_sym)
    
    figure
    hold on
    
    grid on
    box on
    
    plot(K,K,'LineWidth',2)
    
    [nrow, ncol] = size(k_estimates);
    for a = 1:ncol
        plot(K,k_estimates(:,a),k_sym{a},'MarkerSize',5)
    end 
    
    legend(k_names)
    
    xlabel('DPP K (m/s)')
    ylabel('Estimated K (m/s)') 
    
    set(gca,'FontSize',14)

end

