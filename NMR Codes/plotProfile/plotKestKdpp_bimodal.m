function plotKestKdpp_bimodal(K,k_avg,isBimodal)
    
    figure
    hold on
    
    grid on
    box on

    [nrow, ncol] = size(k_avg);
    count = 0;
    
    for a = 1:ncol
        for b = 1:nrow
        %plot(K,k_avg(:,a),'^','MarkerSize',8)
            if isBimodal(b,a) == 1 
%                 scatter(K(b,a),k_avg(b,a),60,'Filled','b','DisplayName','Bimodal')
                count = count + 1;
            else
%                 scatter(K(b,a),k_avg(b,a),30,'r','LineWidth',1.0,'DisplayName','Monomodal')
            end
        end
    end
    
    scatter(K(isBimodal == 1), k_avg(isBimodal == 1),40,'Filled','DisplayName','Bimodal'), 
    scatter(K(isBimodal == 0), k_avg(isBimodal == 0),50,'LineWidth',1.0,'DisplayName','Monomodal'), 

   % legend('Bimodal','Monomodal','Location','northwest')
    legend 
   
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

