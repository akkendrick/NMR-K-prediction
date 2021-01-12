% plot T2 decays
clear
close all

wisc_sites = {'Site1-WellG5','Site1-WellG6','Site2-WellPN1','Site2-WellPN2'};

for k=1:length(wisc_sites)
    site = wisc_sites(k);
   
    [decayCurves,decayTime] = loadRawDecays(site);
    [T2dist, T2logbins, siteName] = loadRawNMRdata(site);
    
    
    % pull out depth intervals
    depths = decayCurves(:,1);
    decayCurves = decayCurves(:,2:end);
    
    decayCurves = fliplr(decayCurves);
    %flipping the data here is it correct?
    %decayTime = fliplr(decayTime);
    
    T2dist = T2dist(:,2:end);


    for k = 1:length(depths)
        depth = depths(k);
        currentT2dist = T2dist(k,:);
        % calculate decay curve from imported data

        E0 = 1;

        T2linbins = 10.^T2logbins; 
        sum = 0;
        for j = 1:length(T2logbins)
            step = currentT2dist(j) .* exp(-decayTime/T2linbins(j));
            sum = sum + step;
            
        end
        
        E_xy{k} = E0 * sum;
        
        % normalize data
        normE_xy{k} = E_xy{k}/max(E_xy{k});
        plotDecay = fliplr(decayCurves(k,:)/max(decayCurves(k,:)));
        
        % make comparison plot
        
        figure(1)
        title('Field decay and inverted decay')
        hold on
        grid on 
        box on
        
        plot(decayTime, plotDecay, 'LineWidth',2)
        plot(decayTime, normE_xy{k},'LineWidth', 2)
        
        xlabel('Time')
        ylabel('Norm Amplitude')
        %set(gca, 'XDir','reverse')
        
        titleString = strcat(string(site),' z= ', string(depths(k)));
        fileString = strcat(string(site),'_z=', string(depths(k)),'.png');
        
        legendStr = ['Measured Decay'; 'Inverted Decay'];
        legend(legendStr,'Location','northeast')
        title(titleString)
        
        print('-dpng','-r300',fileString)
        
        close(1)

    end

end