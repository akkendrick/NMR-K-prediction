%plot all T2 slices 
sites = {'Site1-WellG5','Site1-WellG6','Site2-WellPN1','Site2-WellPN2'};

waterTable = [1.181,1.181,5.0285,4.7476];

cutoffs = [102.0640, 352.0714, 189.5622, 173.5163]; %in ms 
cutoffs = cutoffs*10^-3; %in s

for kk = 1:length(sites)
     baseName = sites{kk}
%      baseDir = '/Volumes/GoogleDrive/My Drive/USGS Project/USGS Data/';
%      gammaBaseDir = '/Volumes/GoogleDrive/My Drive/USGS Project/USGS Data/WI_gamma-EMI-bLS_csvfiles';

    baseDir = 'I:\My Drive\Stanford\USGS Project\Field Data\USGS Data\';
    gammaBaseDir = 'I:\My Drive\Stanford\USGS Project\Field Data\USGS Data\WI_gamma-EMI-bLS_csvfiles';
    
    [T2dist, T2logbins, nmrName] = loadRawNMRdata(sites{kk});
    
   % gammaEMIdepth = [];
   % [gammaEMIdepth, gamma, EMI] = loadGammaEMIData(sites{kk});
    
    [d, K, T2ML, phi, z, SumEch, logK, logT2ML, logPhi, SumEch_3s, SumEch_twm, ...
    SumEch_twm_3s] = loadnmrdata2(nmrName);

    T2depths = T2dist(:,1);
    T2dists{kk} = T2dist(:,2:end);
        
    %%
    zNMR{kk} = T2depths
    
    
end

    T2linbins = 10.^T2logbins;
    
% Organize data for comparision
% zNMR_filt{1} = zNMR{1}(39:63);
% zNMR_filt{2} = zNMR{2}(40:64);

% T2dists_filt{1} = T2dists{1}(39:63,:);
% T2dists_filt{2} = T2dists{2}(40:64,:);

for kk = 1:length(sites)
    for jj = 1:length(zNMR{kk})

        %currentT2dist = T2dists_filt{kk}(:,2:end);
        %currentZ = zNMR_filt{kk};
        
        currentT2dist = T2dists{kk};
        currentZ = zNMR{kk};
        
        T2ML_profile{kk}(jj) = sum(currentT2dist(jj,1:end).*T2logbins)./sum(currentT2dist(jj,1:end));
        
        [val, index] = min(abs(T2ML_profile{kk}(jj) - T2logbins));
        [val2, index2] = min(abs(cutoffs(kk) - T2linbins));
        [val3, index3] = min(abs(33*10^(-3) - T2linbins)); %Compare to 33 ms cutoff

        stepT2dist = currentT2dist(jj,1:end);
        figure(1)
        
        plot(T2linbins, stepT2dist, 'LineWidth',2)
        hold on
        plot(T2linbins(index),currentT2dist(jj,index), '*','MarkerSize',10)
        %plot([T2linbins(index2), T2linbins(index2)],[0, max(stepT2dist(:, index2))], 'r','LineWidth',2)
        %plot([T2linbins(index2), T2linbins(index2)],[0, 0.002], 'm','LineWidth',2)
        plot([T2linbins(index3),T2linbins(index3)], [0,max(stepT2dist)], 'r','LineWidth',2)

        %plot(T2ML_profile{kk}(jj),currentT2dist(jj,index), '*')

        grid on
        box on
        set(gca,'XScale','log')

        hold off
        
        legendString{kk} = strcat('z= ', string(currentZ(jj)),' ',sites{kk});
        fileString = strcat('z= ', string(currentZ(jj)),' ',sites{kk},'.png');
        
        legend(legendString{kk},'Location','northwest')
        print('-dpng','-r300',fileString)
        
    end
    
    
end