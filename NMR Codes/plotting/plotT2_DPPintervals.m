%plot all T2 slices 
sites = {'Site1-WellG5','Site1-WellG6','Site2-WellPN1','Site2-WellPN2'};

waterTable = [1.181,1.181,5.0285,4.7476];
depthOffsets = [0.95,0.75,0.75,0.75];

% 
% cutoffs = [102.0640, 352.0714, 189.5622, 173.5163]; %in ms (median k diff factor)
% cutoffs = cutoffs*10^-3; %in s

for kk = 1:length(sites)
     baseName = sites{kk}
%      baseDir = '/Volumes/GoogleDrive/My Drive/USGS Project/USGS Data/';
%      gammaBaseDir = '/Volumes/GoogleDrive/My Drive/USGS Project/USGS Data/WI_gamma-EMI-bLS_csvfiles';

    baseDir = 'I:\My Drive\Stanford\USGS Project\Field Data\USGS Data\';
    gammaBaseDir = 'I:\My Drive\Stanford\USGS Project\Field Data\USGS Data\WI_gamma-EMI-bLS_csvfiles';

    if strcmp(baseName,'Site1-WellG5')
        site = 'Site1-WellG5';
        name = 'G5_W1_tr5_20x_16p5_up_F1n2_wRIN_wRFI_Reg50_Va1';
        nmrName = name;
    elseif strcmp(baseName,'Site1-WellG6')
        site = 'Site1-WellG6';
        name = 'G6_W2_tr5_20x_16p75_up_F_wRIN_wRFI_reg50_Va1';
        nmrName = name;
    elseif strcmp(baseName,'Site2-WellPN1')
        site = 'Site2-WellPN1';
        name = 'Pl_W1_Tr5_20x_MPp75aLS_F1n2_wRIN_wRFI_Reg50_Va1';
        nmrName = name;
    elseif strcmp(baseName,'Site2-WellPN2')
        site = 'Site2-WellPN2';
        name = 'W2_Tr5_20x_MPp75aLS_Reg50_wRIN_wRFI_Va1';
        nmrName = name;
    end
    
    in1 = [baseDir site '/' site '_SE_decay' '.txt']; 
    in2 = [baseDir site '/' site '_1Dvectors.txt'];
    in3 = [baseDir site '/' strcat(site,'_DPP_filt.txt')];
    in4 = [baseDir site '/' site '_SE_decay_time.txt'];
    in5 = [baseDir site '/' site '_T2_dist.txt'];
    in6 = [baseDir site '/' site  '_T2_bins_log10s.txt'];
    
    decaycurve = load(in1); 
    dparam = dlmread(in2,'\t',1,0); 
    DPPdat = load(in3); 
    t = load(in4);
    
    Dk = DPPdat(:,2)*1.16e-5; % converts K from m/day to m/s
    z_dk = DPPdat(:,1);
    
    T2dist = load(in5);
    T2dist(:,1) = T2dist(:,1) - depthOffsets(kk);
    
    T2dists{kk} = T2dist;
    
    T2depths = T2dist(:,1);
    T2logbins = load(in6);
    
    % set needed variables
    NMRphi = dparam(:,2);
    S = decaycurve(:,2:end); 

    [d, K, T2ML, phi, z, SumEch, logK, logT2ML, logPhi, SumEch_3s, SumEch_twm, ...
    SumEch_twm_3s] = loadnmrdata2(baseName); 

    origT2dist = T2dist;
    
    errorz = zeros(1,length(Dk));
    ind = zeros(1,length(Dk));
    
    for i = 1:length(Dk)
       [errorz(i), ind(i)] = min(abs(z_dk(i) - T2depths)); 
    end


    %%
    zNMR = T2depths(ind);
    
    T2linbins = 10.^T2logbins;
    
    for jj = 1:length(zNMR)
        
        %currentT2dist = T2dists_filt{kk}(:,2:end);
        %currentZ = zNMR_filt{kk};
        T2distIndex = ind(jj);
        
        currentT2dist = T2dists{kk};
        currentZ = zNMR;
        
        T2ML_profile{kk}(jj) = sum(currentT2dist(T2distIndex,2:end).*T2logbins)./sum(currentT2dist(T2distIndex,2:end));
        
        [val, index] = min(abs(T2ML_profile{kk}(jj) - T2logbins));
%         [val2, index2] = min(abs(cutoffs(kk) - T2linbins));
        [val3, index3] = min(abs(33*10^(-3) - T2linbins)); %Compare to 33 ms cutoff

        stepT2dist = currentT2dist(T2distIndex,2:end);
        figure(1)
        
        plot(T2linbins, stepT2dist, 'LineWidth',2)
        hold on
        plot(T2linbins(index),currentT2dist(T2distIndex,index), '*','MarkerSize',10)
        %plot([T2linbins(index2), T2linbins(index2)],[0, max(stepT2dist(:, index2))], 'r','LineWidth',2)
%         plot([T2linbins(index2), T2linbins(index2)],[0, 0.002], 'm','LineWidth',2)
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