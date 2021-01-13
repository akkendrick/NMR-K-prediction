% plot T2 decays
clear
close all

wisc_sites = {'Site1-WellG5','Site1-WellG6','Site2-WellPN1','Site2-WellPN2'};

for m=1:length(wisc_sites)
    site = wisc_sites(m);
   
    [decayCurves,decayTime] = loadRawDecays(site);
    [T2dist, T2logbins, siteName] = loadRawNMRdata(site);
    
    [d, K, T2ML, phi{m}, z, SumEch, logK, logT2ML, logPhi, SumEch_3s, SumEch_twm, ...
    SumEch_twm_3s] = loadnmrdata2(site{1});     
    
    % pull out depth intervals
    depths = d;
    decayCurves = decayCurves(:,2:end);
    
    decayCurves = fliplr(decayCurves);
    %flipping the data here is it correct?
    %decayTime = fliplr(decayTime);
    
    T2dist = T2dist(:,2:end);
    sitePhi = phi{m};


    for k = 1:length(depths)
        depth = depths(k);
        currentT2dist = T2dist(k,:);
        % calculate decay curve from imported data


        T2linbins = 10.^T2logbins; 
        sumInv = 0;
        for j = 1:length(T2logbins)
            step = currentT2dist(j) .* exp(-decayTime/T2linbins(j));
            sumInv = sumInv + step;
            
        end
        
        plotDecay = fliplr(decayCurves(k,:));
        
        sumInv = sumInv./max(sumInv);
        
        E0 = sitePhi(k)*100;
        E_xy{k} = E0 * sumInv;
        
        % normalize data
        %normE_xy{k} = E_xy{k}/max(E_xy{k});
        
        
        
        % make comparison plot
        
        figure(1)
        title('Field decay and inverted decay')
        hold on
        grid on 
        box on
        
        plot(decayTime, plotDecay, 'LineWidth',2)
        plot(decayTime, E_xy{k},'LineWidth', 2)
        
        xlabel('Time')
        ylabel('Amplitude')
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