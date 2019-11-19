function plotKPhiwithDepth(K,NMRphi,waterTable,z,T2dist,T2logbins,k_estimates,k_names,k_sym,noPorosity,bProfile)
% Plot 


if noPorosity == 1
figure('Renderer', 'painters', 'Position', [10 10 600 900])

    %T2dist = flip(T2dist);

    depths = T2dist(:,1);
    T2dist = T2dist(:,2:end);
      
    subplot(121)
    %imagesc(T2dist)
    hold on
    
    waterTableLine = ones(1,length(T2logbins));
    waterTableLine = waterTableLine .* waterTable;
    
    surf(T2logbins,depths,T2dist,'EdgeColor','none')
    plot3(T2logbins,waterTableLine,ones(1,length(T2logbins)),'r','LineWidth',3)
    
    view(0,90)
    set(gca, 'YDir','reverse')
    ylim([min(depths),max(depths)])
    xlim([min(T2logbins),max(T2logbins)])
    ylabel('Depth (m)')
    xlabel('log_{10} T_2 (s)')
    set(gca,'FontSize',16)
  %   set(gca,'xscale','log')
    
  %  set(gca, 'YDir','reverse')
  
    
    subplot(122)
    
    waterTableLine = [waterTable waterTable];
    
    box on
    grid on
    
    hold on
    plot(K,z,'*','MarkerSize',8)

    [nrow, ncol] = size(k_estimates);
    for a = 1:ncol
        plot(k_estimates(:,a),z,k_sym{a},'MarkerSize',8)
    end 
       
    plot([0,1],waterTableLine,'r','LineWidth',3,'HandleVisibility','off')
    
    legend([{'DPP'} k_names])

   
    xlabel(strcat('\it K', '\rm (m/s)'))
    set(gca, 'YDir','reverse')
    ylim([min(depths),max(depths)])
    %xlim([min(K),max(K)])

    xlim([10^-6, 10^-3])
    set(gca,'XScale','log')
    
    set(gca,'FontSize',16)
    
else
    
    figure

    %T2dist = flip(T2dist);

    depths = T2dist(:,1);
    T2dist = T2dist(:,2:end);
      
    subplot(141)
    %imagesc(T2dist)
    hold on
    
    waterTableLine = ones(1,length(T2logbins));
    waterTableLine = waterTableLine .* waterTable;
    
    surf(T2logbins,depths,T2dist,'EdgeColor','none')
    plot3(T2logbins,waterTableLine,ones(1,length(T2logbins)),'r','LineWidth',3)
    
    view(0,90)
    set(gca, 'YDir','reverse')
    ylim([min(depths),max(depths)])
    xlim([min(T2logbins),max(T2logbins)])
    ylabel('Depth (m)')
    xlabel('log_{10} T_2 (s)')
    set(gca,'FontSize',16)
  %   set(gca,'xscale','log')
    
  %  set(gca, 'YDir','reverse')
  
    
    subplot(142)
    
    waterTableLine = [waterTable waterTable];
       
    box on
    grid on
    
    hold on
    plot(K,z,'*','MarkerSize',8)

    [nrow, ncol] = size(k_estimates);
    for a = 1:ncol
        plot(k_estimates(:,a),z,k_sym{a},'MarkerSize',5)
    end 
    
    %legend(k_names)
    
    plot([10^-6,10^-2],waterTableLine,'r','LineWidth',3,'HandleVisibility','off')

   
    xlabel(strcat('\it K', '\rm (m/s)'))
    set(gca, 'YDir','reverse')
    ylim([min(depths),max(depths)])
    set(gca,'XScale','log')
    set(gca,'FontSize',16)

    xlim([10^-6, 10^-2])

    subplot(143)
    
    waterTableLine = [waterTable waterTable];
    
    box on
    grid on
    
    hold on
    plot(NMRphi,depths,'LineWidth',3)
    plot([0,1],waterTableLine,'r','LineWidth',3)

    
    xlabel('NMR Porosity')
    set(gca, 'YDir','reverse')
    ylim([min(depths),max(depths)])
    xlim([0,0.5])
    
    set(gca,'FontSize',16)
    
    subplot(144)
    grid on
    box on
    hold on
    waterTableLine = [waterTable waterTable];

    
    plot([0,1],waterTableLine,'r','LineWidth',3,'HandleVisibility','off')
    scatter(bProfile,z,30,'Filled')
    set(gca,'YDir','reverse')
    xlim([min(bProfile),max(bProfile)])
    ylim([min(depths),max(depths)])
    
end


end

