function plotKPhiwithDepth_v2(K,NMRphi,waterTable,z,T2dist,T2logbins,k_est,k_avg,k_names,noPorosity,plotIndex)
% Plot 

[z, index] = sort(z,'descend');

K = K(index);
k_est = k_est(index,:);
k_avg = k_avg(index,:);

minDepth = 1;
maxDepth = 18;

f1 = figure('Renderer', 'painters', 'Position', [10 10 300 300]);
figure(f1)


if noPorosity == 1
    

    %T2dist = flip(T2dist);

    depths = T2dist(:,1);
    T2dist = T2dist(:,2:end);
      
    subplot(1,2,1)
    %imagesc(T2dist)
    hold on
    
    waterTableLine = ones(1,length(T2logbins));
    waterTableLine = waterTableLine .* waterTable;
    
    surf(T2logbins,depths,T2dist,'EdgeColor','none')
    plot3(T2logbins,waterTableLine,ones(1,length(T2logbins)),'m','LineWidth',3)
    
    view(0,90)
    set(gca, 'YDir','reverse')
    ylim([minDepth,maxDepth])
    xlim([min(T2logbins),max(T2logbins)])
    ylabel('Depth (m)')
    xlabel('log_{10} T_2 (s)')
    set(gca,'FontSize',12)
  %   set(gca,'xscale','log')
    grid on
    box on
    hold off
  
  %  set(gca, 'YDir','reverse')
  
    
    subplot(1,2,2)
    
    waterTableLine = [waterTable waterTable];
    
    box on
    grid on
    
    hold on

    color1 = [222 235 247]/255;
    color2 = [49 130 189]/255;
    numcolors = length(z);
    
    % create the gradients
    clrmap = cell2mat(arrayfun(@(a,b)linspace(a,b,numcolors)',color2,color1,'uni',false));
    
    [nrow, ncol] = size(k_avg)
    for a = 1:ncol
        %plot(k_avg(:,a),z,'^','MarkerSize',8)
        scatter(k_avg(:,a),z,40,'filled')
    end 
    plot(K,z,'or','MarkerSize',8)

    %legend([{'DPP'} k_names])
 
%     [nrow, ncol] = size(k_est);
%     for a = 1:ncol
%         plot(k_est(:,a),z,'k.','MarkerSize',2,'HandleVisibility','off')
%     end 
        
    plot([0,1],waterTableLine,'m','LineWidth',3,'HandleVisibility','off')

   
    xlabel(strcat('\it K', '\rm (m/s)'))
    set(gca, 'YDir','reverse')
    ylim([minDepth,maxDepth])
    %xlim([min(K),max(K)])
%     if plotIndex > 4
%         xlim([10^-5, 10^-3])
%     else
%         xlim([10^-6, 10^-3])
%     end

    xlim([10^-6, 10^-2])
    set(gca,'XScale','log')
    
    set(gca,'FontSize',12)
    
    scale = 0.05;
    pos = get(gca, 'Position');
    pos(2) = pos(2)+scale*pos(4);
    pos(4) = (1-scale)*pos(4);
    set(gca, 'Position', pos)
    hold off
else
    
    figure

        grid on
    box on
    %T2dist = flip(T2dist);

    depths = T2dist(:,1);
    T2dist = T2dist(:,2:end);
      
    subplot(131)
    %imagesc(T2dist)
    hold on
    grid on
    box on
    waterTableLine = ones(1,length(T2logbins));
    waterTableLine = waterTableLine .* waterTable;
    
    surf(T2logbins,depths,T2dist,'EdgeColor','none')
    plot3(T2logbins,waterTableLine,ones(1,length(T2logbins)),'m','LineWidth',3)
    
    view(0,90)
    set(gca, 'YDir','reverse')
    ylim([minDepth,maxDepth])
    xlim([min(T2logbins),max(T2logbins)])
    ylabel('Depth (m)')
    xlabel('log_{10} T_2 (s)')
    set(gca,'FontSize',16)
    grid on
    box on
  %   set(gca,'xscale','log')
    
  %  set(gca, 'YDir','reverse')
  
    
    subplot(132)
    
    waterTableLine = [waterTable waterTable];
    
    box on
    grid on
    
    hold on
    plot(K,z,'*','MarkerSize',8)

    [nrow, ncol] = size(k_est);
    for a = 1:ncol
        plot(k_est(:,a),z,k_sym{a},'MarkerSize',5)
    end 
    
   legend(k_names)
    
   plot([0,1],waterTableLine,'m','LineWidth',3)

   
    xlabel(strcat('\it K', '\rm (m/s)'))
    set(gca, 'YDir','reverse')
    ylim([minDepth,maxDepth])
    xlim([min(K),max(K)])
    set(gca,'FontSize',16)

    
    
    subplot(133)
    
    waterTableLine = [waterTable waterTable];
    
    box on
    grid on
    
    hold on
    plot(NMRphi,depths,'LineWidth',3)
    plot([0,1],waterTableLine,'m','LineWidth',3)

    
    xlabel('NMR Porosity')
    set(gca, 'YDir','reverse')
    ylim([minDepth,maxDepth])
    xlim([0,1])
    
    set(gca,'FontSize',16)
    
end


end

