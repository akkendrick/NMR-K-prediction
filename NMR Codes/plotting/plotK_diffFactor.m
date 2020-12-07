function plotK_diffFactor(K,k_estimates,k_avg,k_names)
    
    %figure
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
   
    magLine = ones(length(k_avg(:,1))).*10;
    plot([10^-7,10^-2],[10,10],'k','HandleVisibility','off','LineWidth',2)
    
    xlabel('DPP K (m/s)')
    ylabel('K Difference Factor') 
    
    set(gca,'FontSize',16)
    set(gca,'XScale','log')
    set(gca,'YScale','log')
    
     %ylim([10^-6,10^-2])
     xlim([10^-7,10^-2])
end

