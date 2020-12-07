function plotKestKdpp_sym(K,k_estimates,k_avg,k_names,k_sym,k_size)
    
    figure
    hold on
    
    grid on
    box on
%     smoothedData = smooth(k_avg,10,'lowess');
    [nrow, ncol] = size(k_avg);
    for a = 1:ncol
        %plot(K,k_avg(:,a),'^','MarkerSize',8)
        plot(K,k_avg(:,a),k_sym(a),'MarkerSize',k_size(a))
        plot(K,smooth(k_avg(:,a),8,'lowess'),'LineWidth',2,'HandleVisibility','off')

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

