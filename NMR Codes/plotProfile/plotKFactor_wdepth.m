function plotKFactor_wdepth(K,k_est,k_avg,z,k_names)
    

    [z, index] = sort(z, 'descend');
    
    K = K(index);
    k_est = k_est(index,:);
    k_avg = k_avg(index,:);


    kDiffFactor = estimateKdiffFactor(K,k_est,0);
    kDiffFactorAvg = estimateKdiffFactor(K,k_avg,0);
    

    
   % kUnitFactor = estimateKdiffFactor(K,K);
    
    hold on
    
    grid on
    box on
    
%     color1 = [222 235 247]/255;
%     color2 = [49 130 189]/255;
    
    color3 = [253,224,221]/255;
    color4 = [197,27,138]/255;


    color1 = [207 218 155]/255;
    color2 = [44 127 184]/255;
    
    numcolors = length(z);

    % create the gradients
    clrmap = cell2mat(arrayfun(@(a,b)linspace(a,b,numcolors)',color2,color1,'uni',false));
    clrmap2 = cell2mat(arrayfun(@(a,b)linspace(a,b,numcolors)',color4,color3,'uni',false));

%     [nrow, ncol] = size(k_avg);
%     for a = 1:ncol
%         scatter(K,kDiffFactorAvg(:,a),60,clrmap,'filled')
%         %scatter(K,kDiffFactorAvg(:,a),40,'filled')
%     end 
    
%     scatter(K,kDiffFactorAvg(:,1),60,clrmap,'filled')
%     scatter(K,kDiffFactorAvg(:,2),60,clrmap2,'filled')
    
    scatter(K,kDiffFactorAvg(:,1),60,'filled')
%    legend(k_names,'Location','northwest')

%     [nrow, ncol] = size(k_est);
%     for a = 1:ncol
%         plot(K,kDiffFactor(:,a),'.k','MarkerSize',2,'HandleVisibility','off')
%     end 
    
    plot(K,ones(1,length(K)),'k','LineWidth',2,'HandleVisibility','off')
  
    
    xlabel('DPP K (m/s)')
    ylabel('K Difference Factor') 
    
    set(gca,'FontSize',14)
    set(gca,'XScale','log')
    set(gca,'YScale','log')
    

end

