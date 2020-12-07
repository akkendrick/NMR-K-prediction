function plotKestKdpp_v2(K,k_estimates,k_avg,k_names)
    
    figure
    hold on
    
    grid on
    box on

    [nrow, ncol] = size(k_avg);
    for a = 1:ncol
        %plot(K,k_avg(:,a),'^','MarkerSize',8)
        scatter(K,k_avg(:,a),60,'Filled')

    end 
    
    legend(k_names,'Location','northwest')
         
    [nrow, ncol] = size(k_estimates);
    for a = 1:ncol
        %plot(K,k_estimates(:,a),'k.','MarkerSize',2,'HandleVisibility','off')
    end 
    
    plot(K,K,'k','LineWidth',2,'HandleVisibility','off')
    plot(K,K*10,'k:','HandleVisibility','off')
    plot(K,K*0.1,'k:','HandleVisibility','off')
    
    xlabel('DPP K (m/s)')
    ylabel('Estimated K (m/s)') 
    
    set(gca,'FontSize',16)
    set(gca,'XScale','log')
    set(gca,'YScale','log')
    
%     ylim([10^-6,10^-2])
%     xlim([4*10^-6,10^-3])
end

