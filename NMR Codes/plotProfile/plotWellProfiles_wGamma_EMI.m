function plotWellProfiles_wGamma_EMI(K,NMRphi,waterTable,z,T2ML,T2dist,T2logbins,...
    gammaEMIdata, gammaEMIoffset, DPPdat)
% Plot 
    
    figure('Renderer', 'painters', 'Position', [10 10 800 500])


    %T2dist = flip(T2dist);
    minDepth = 1;
    maxDepth = 18;
    depths = T2dist(:,1);
    T2dist = T2dist(:,2:end);
      
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot T2 distribution
    
    % Compute T2ML for all points
    for kk = 1:length(depths)
        T2ML_profile(kk) = sum(T2dist(kk,:).*T2logbins)./sum(T2dist(kk,:));
    end
    
    subplot(143)
    %imagesc(T2dist)
    hold on
    box on
    grid on
    
    waterTableLine = ones(1,length(T2logbins));
    waterTableLine = waterTableLine .* waterTable;
    
    surf(T2logbins,depths,T2dist,'EdgeColor','none')
    %plot3(T2logbins,waterTableLine,ones(1,length(T2logbins)),'-.m','LineWidth',3)
    %plot3(log10(T2ML), z, ones(1,length(T2ML)),'*r','MarkerSize',8)
    plot3(T2ML_profile,depths,ones(1,length(T2ML_profile)),'*r','MarkerSize',8)
%     for kk = 1:length(depths)
%         depths(kk)
%         currentDepth = ones(length(T2dist),1)*depths(kk);
%         plot3(T2logbins',currentDepth,T2dist)
%     
%     end
        
    view(0,90)
    set(gca, 'YDir','reverse')
    %ylim([minDepth,maxDepth])
    ylim([minDepth,maxDepth])

    xlim([min(T2logbins),max(T2logbins)])
    set(gca,'YTickLabel',[]);
    xlabel('NMR log_{10} T_2 (s)')

    set(gca,'FontSize',14)
  %   set(gca,'xscale','log')
    
  %  set(gca, 'YDir','reverse')

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
   % Plot lithologies 
   
    handle = subplot(141);
    set(gca,'YTickLabel',[]);
    set(gca,'XTickLabel',[]);
    set(gca,'XTick',[]);
    set(gca,'YTick',[]);
      
    origPosition = get(handle,'Position');
   
    gammaEMIdepth = gammaEMIdata(:,1) - gammaEMIoffset;
    gamma = gammaEMIdata(:,2);
    EMI = gammaEMIdata(:,3);

    %We have already corrected for depth here! See above 
    EMIdepth = gammaEMIdepth(EMI > 0);
    goodEMI = EMI(EMI > 0);
   
%     origPosition = get(0,'defaultfigureposition');
%     ax1 = axes('Position',origPosition);
     ax1 = axes('Position',origPosition,'Color','none');
     hold(ax1, 'on')

    smoothWindow = 10;
    smoothGamma = smooth(gamma, smoothWindow);

    %plot(gamma, depth, 'parent', ax1)
    plot(smoothGamma, gammaEMIdepth,'parent',ax1,'LineWidth',2)

    ylim([minDepth,maxDepth])
    xlim([0,100])
    
    set(ax1,'YDir','reverse')

    ax1.XColor = 'k';
    ax1.YColor = 'k';

    xlabel('Gamma (Counts/second)')
    ylabel('Depth (m)')
    set(gca,'FontSize',14)

    box on 
    grid on
    grid minor

    ax1_pos = ax1.Position; % position of first axes
    ax2 = axes('Position',ax1_pos,...
        'XAxisLocation','top',...
        'YAxisLocation','right',...
        'Color','none','NextPlot','add');
    ax2.XColor = 'r';

    hold(ax2, 'on')


    plot(goodEMI, EMIdepth, '--r', 'parent', ax2,'LineWidth',2)
    set(gca,'YDir','reverse')
    set(gca,'YTickLabel',[])
    
    ylim([minDepth,maxDepth])
    xlim([0,100])
    set(gca,'FontSize',14)

    xlabel('EMI Conductivity (mS/m)')

  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    % Plot K Estimates
    handle = subplot(144);
    set(gca,'YTickLabel',[]);
    set(gca,'XTickLabel',[]);
    set(gca,'XTick',[]);
    set(gca,'YTick',[]);
    
    origPosition = get(handle,'Position');

    DPP_depths = DPPdat(:,1);
    DPP_K = DPPdat(:,2)*1.12e-5;
    
    waterTableLine = [waterTable waterTable];
    hold on
    
    %scatter(K,z,60,'Filled')

%     [nrow, ncol] = size(k_estimates);
%     for a = 1:ncol
%         plot(k_estimates(:,a),z,k_sym{a},'MarkerSize',5)
%     end
    
    %scatter(dlubacModel,z,40,[0.9290,0.6940,0.1250],'Filled')
    %scatter(SDRModel,z,40,[0.466,0.6740,0.1880],'Filled')
    
    ax1 = axes('Position',origPosition,'Color','none');
    hold(ax1, 'on')
    
    plot(DPP_K, DPP_depths,'or','parent',ax1,'MarkerSize',8)
    
    box on
    grid on
    
    %plot(K,z,'og','MarkerSize',10)

    %legend(k_names)
    
    %plot([10^-6,10^-2],waterTableLine,'r','LineWidth',3,'HandleVisibility','off')

   
    xlabel(strcat('DPP \it K', '\rm (m/s)'))
    set(gca, 'YDir','reverse')
    ylim([minDepth,maxDepth])
   
    set(gca,'XScale','log')
%    set(gca,'FontSize',14)
%    set(gca,'YTickLabel',[]);
% 
%     xh = get(gca,'xlabel');
%     p = get(xh,'position'); % get the current position property
%     p(2) = 2*p(2) ;        % double the distance, 
%                            % negative values put the label below the axis
%     set(xh,'position',p) 
        set(gca,'FontSize',14)

    xlim([10^-6, 10^-2])
    
%     xMin_md = 10^-6 * 86400;
%     xMax_md = 10^-2 * 86400;
%     
%     ax1_pos = ax1.Position; % position of first axes
%     ax2 = axes('Position',ax1_pos,...
%         'XAxisLocation','top',...
%         'YAxisLocation','right',...
%         'Color','none','NextPlot','add');
%     ax2.XColor = 'r';

   % hold(ax2, 'on')
    
    %DPP_K_mday = DPP_K * 86400;
    
    %plot(DPP_K_mday, DPP_depths, 'or', 'parent', ax2,'MarkerSize',8)
    
    %set(gca,'YDir','reverse')
    %set(gca,'YTickLabel',[])
    
    %ylim([minDepth,maxDepth])
    %xlim([xMin_md,xMax_md])
    
    %set(gca,'FontSize',12)
    %set(gca, 'XScale','log')
    %xlabel(strcat('\it K', '\rm (m/d)'))


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot Porosity
    subplot(142)
    waterTableLine = [waterTable waterTable];
    
    box on
    grid on
    grid minor
    
    hold on
    plot(NMRphi,depths,'LineWidth',3)
    %plot([0,1],waterTableLine,'r','LineWidth',3)
     
    
    xlabel('NMR WC')
    set(gca, 'YDir','reverse')
    ylim([minDepth,maxDepth])
    xlim([0,0.6])
%    set(gca,'FontSize',14)
   set(gca,'YTickLabel',[]);
   set(gca, 'XMinorTick', 'on')
    set(gca,'FontSize',14)

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
%     ylim([minDepth,maxDepth])
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
%     ylim([minDepth,maxDepth])
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
%     ylim([minDepth,maxDepth])
%     xlim([10^-4 10^0])
%     set(gca,'XScale','log')
%     set(gca,'YTickLabel',[]);
%     xlabel('SDR b')
%     set(gca,'FontSize',12)

    
end




