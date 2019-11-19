function plotKestKdpp_wdepth(K,k_est,k_avg,z,k_names,colors)
    
    [z, index] = sort(z,'descend');

    K = K(index);
    k_est = k_est(index,:);
    k_avg = k_avg(index,:);

    
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

    [nrow, ncol] = size(k_avg);
    for a = 1:ncol
        %scatter(K,k_avg(:,a),60,clrmap,'filled')
        %scatter(K,k_avg(:,a),40,colors{a},'filled')
        scatter(K,k_avg(:,a),40,'filled')
    end 

    %scatter(K,k_avg(:,1),60,clrmap,'filled')
    %scatter(K,k_avg(:,2),60,clrmap2,'filled')
    
%     scatter(K,k_avg(:,1),60,'filled')
    
    
%    legend(k_names,'Location','northwest')
         
    [nrow, ncol] = size(k_est);
    for a = 1:ncol
       % plot(K,k_est(:,a),'k.','MarkerSize',2,'HandleVisibility','off')
    end 
    
    plot(K,K,'k','LineWidth',2,'HandleVisibility','off')
    plot(K,K*10,'k:','HandleVisibility','off')
    plot(K,K*0.1,'k:','HandleVisibility','off')
    
    xlabel('DPP K (m/s)')
    ylabel('Estimated K (m/s)') 
    
    set(gca,'XScale','log')
    set(gca,'YScale','log')
    set(gca,'FontSize',12)
    
    ylim([10^-6,10^-2])
    %xlim([4*10^-6,10^-3])
    xlim([4*10^-6,10^-3])
end

