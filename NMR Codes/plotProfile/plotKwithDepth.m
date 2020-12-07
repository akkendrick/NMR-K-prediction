function plotKwithDepth(K,z,T2dist,T2logbins,k_estimates,k_names,k_sym)
% Plot 

    
    figure

    %T2dist = flip(T2dist);

    depths = T2dist(:,1);
    T2dist = T2dist(:,2:end);
      
    subplot(121)
    %imagesc(T2dist)
    
    surf(T2logbins,depths,T2dist,'EdgeColor','none')
    
    view(0,90)
    set(gca, 'YDir','reverse')
    ylim([min(depths),max(depths)])
    xlim([min(T2logbins),max(T2logbins)])
    ylabel('Depth (m)')
    xlabel('log_{10} T_2 (s)')
    set(gca,'FontSize',14)
  %   set(gca,'xscale','log')
    
  %  set(gca, 'YDir','reverse')
  
    
    subplot(122)
    
    box on
    grid on
    
    hold on
    plot(K,z,'*','MarkerSize',8)

    [nrow, ncol] = size(k_estimates);
    for a = 1:ncol
        plot(k_estimates(:,a),z,k_sym{a},'MarkerSize',5)
    end 
    
    legend(k_names)
    
    xlabel('Hydraulic Conductivity (m/s)')
    set(gca, 'YDir','reverse')
    ylim([min(depths),max(depths)])
    set(gca,'XScale','log')
    set(gca,'FontSize',14)

      
    
end

