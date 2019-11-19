function plotWellProfiles_wLith(K,NMRphi,waterTable,z,T2ML,T2dist,T2logbins,...
    k_estimates,k_names,k_sym,bProfile,lithStart,lithEnd,lithUnits,SDRModel,dlubacModel,meanT2B)
% Plot 
    
    figure('Renderer', 'painters', 'Position', [10 10 600 900])


    %T2dist = flip(T2dist);

    depths = T2dist(:,1);
    T2dist = T2dist(:,2:end);
      
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot T2 distribution
    
    subplot(142)
    %imagesc(T2dist)
    hold on
    
    waterTableLine = ones(1,length(T2logbins));
    waterTableLine = waterTableLine .* waterTable;
    
    surf(T2logbins,depths,T2dist,'EdgeColor','none')
    plot3(T2logbins,waterTableLine,ones(1,length(T2logbins)),'r','LineWidth',3)
    plot3(log10(T2ML), z, ones(1,length(T2ML)),'*r','MarkerSize',8)
    
    view(0,90)
    set(gca, 'YDir','reverse')
    ylim([min(depths),max(depths)])
    xlim([min(T2logbins),max(T2logbins)])
    set(gca,'YTickLabel',[]);
    xlabel('log_{10} T_2 (s)')
   % set(gca,'FontSize',14)
  %   set(gca,'xscale','log')
    
  %  set(gca, 'YDir','reverse')

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
   % Plot lithologies 
   subplot(141)
   hold on
   grid on
   box on
   
    clay = [128,128,0]/255;
    sand = [218,165,32]/255;
    gravel = [169,169,169]/255;
    bedrock = [0,0,0]/255;
    soil = [100,125,150]/255;
    
    colors = [clay; sand; gravel; bedrock; soil];
    
    startPosition = 1;
    width = 5;

    for i=1:length(lithStart)
        if lithStart(i) == 0
           startPosition = 0;
           rectangle('Position',[startPosition,-1*lithEnd(i),width,abs(lithEnd(i))-lithStart(i)],'FaceColor',colors(lithUnits(i),:))
        else
           rectangle('Position',[startPosition,-1*lithEnd(i),width,abs(lithEnd(i))-lithStart(i)],'FaceColor',colors(lithUnits(i),:))
        end

    end

    xlabel('Known Lithology')
    %set(gca,'YDir','reverse')

%   set(gca,'FontSize',14)
    ylim([-max(depths),-min(depths)])
    ylabel('Depth (m)')
   set(gca,'XTickLabel',[])
    
    set(gca, 'YTickLabel',[16 14 12 10 8 6 4 2])
%     subplot(172)
%     
%     waterTableLine = [waterTable waterTable];
%     
%     box on
%     grid on
%     hold on
%     plot([-1,0],waterTableLine,'r','LineWidth',3,'HandleVisibility','off')
% 
%     scatter(log10(T2ML),z,30,'Filled')
%     
%     
%     set(gca,'FontSize',16)
%     set(gca, 'YDir','reverse')
%     ylim([min(depths),max(depths)])
%     xlabel('log_{10} T_{2ML}')
%   
  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    % Plot K Estimates
    subplot(144)
    
    waterTableLine = [waterTable waterTable];
    box on
    grid on
    hold on
    
    %scatter(K,z,60,'Filled')

%     [nrow, ncol] = size(k_estimates);
%     for a = 1:ncol
%         plot(k_estimates(:,a),z,k_sym{a},'MarkerSize',5)
%     end
    
    %scatter(dlubacModel,z,40,[0.9290,0.6940,0.1250],'Filled')
    %scatter(SDRModel,z,40,[0.466,0.6740,0.1880],'Filled')
    plot(K,z,'or','MarkerSize',10)


    %legend(k_names)
    
    %plot([10^-6,10^-2],waterTableLine,'r','LineWidth',3,'HandleVisibility','off')

   
    xlabel(strcat('\it K', '\rm (m/s)'))
    set(gca, 'YDir','reverse')
    ylim([min(depths),max(depths)])
   
    set(gca,'XScale','log')
%    set(gca,'FontSize',14)
%    set(gca,'YTickLabel',[]);
% 
%     xh = get(gca,'xlabel');
%     p = get(xh,'position'); % get the current position property
%     p(2) = 2*p(2) ;        % double the distance, 
%                            % negative values put the label below the axis
%     set(xh,'position',p) 
    
    xlim([10^-6, 10^-2])
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot Porosity
    subplot(143)
    waterTableLine = [waterTable waterTable];
    
    box on
    grid on
    
    hold on
    plot(NMRphi,depths,'LineWidth',3)
    %plot([0,1],waterTableLine,'r','LineWidth',3)
     
    
    xlabel('NMR Porosity')
    set(gca, 'YDir','reverse')
    ylim([min(depths),max(depths)])
    xlim([0,0.6])
%    set(gca,'FontSize',14)
   set(gca,'YTickLabel',[]);

%     ax1 = gca;
%     ax1.XColor = 'r';
%     ax1_pos = ax1.Position;
%         
%     ax2 = axes('Position',ax1_pos,'XAxisLocation','top',...
%         'YAxisLocation','right','Color','none');
%     hold on
%     plot(log10(T2ML), z, '*','Parent',ax2)
%     
%     box on; grid on
%     set(gca,'YDir','reverse')
%     ylim([min(depths),max(depths)])
%     xlim([-3, 0])
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot T2ml
%     subplot(154)
%     grid on
%     box on
%     hold on
%     
%     T2Bcutoff = meanT2B*0.33333; % Dlubac 2014 says T2ML > 1/3 T2B indicates we should care about T2B
%     T2Bprofile = ones(1,length(z)).*T2Bcutoff;
%     
%     scatter(log10(T2ML),z,40,'filled')
%   % scatter(log10(T2Bprofile),z,40,'filled')
%     
%     set(gca,'YDir','reverse')
%     ylim([min(depths),max(depths)])
%     xlim([-3, 0])
%     xlabel('log_{10} T_{2ML}')
% %    set(gca,'FontSize',14)
%    set(gca,'YTickLabel',[]);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % Plot b values
%     subplot(154)
%     grid on
%     box on
%     hold on
%     waterTableLine = [waterTable waterTable];
% 
%     
%     %plot([0,1],waterTableLine,'r','LineWidth',3,'HandleVisibility','off')
%     scatter((bProfile),z,30,'Filled')
%     
%     set(gca,'YDir','reverse')
%     %xlim([min(bProfile),max(bProfile)])
%     ylim([min(depths),max(depths)])
%     xlim([10^-4 10^0])
%     set(gca,'XScale','log')
%     set(gca,'YTickLabel',[]);
%     xlabel('SDR b')
%     set(gca,'FontSize',16)

    
end




