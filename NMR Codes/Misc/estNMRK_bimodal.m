% Estimate Krange using bimodal data

clear 

%sites = {'Site1-WellG5','Site1-WellG6','Site2-WellPN1','Site2-WellPN2'};

sites = {'Site1-WellG6'};
legendNames = [{'Site1-WellG6'}];

waterTable = [2.1248,2.0469,5.0285,4.7476]; % rel ground surface
depthOffsets = [0.95,0.75,0.75,0.75];

for kk = 1:length(sites)
    %close all
    baseName = sites{kk}
%     
    rawBaseDir = '/Volumes/GoogleDrive/My Drive/USGS Project/USGS Data/';
    baseDir = '/Volumes/GoogleDrive/My Drive/USGS Project/NMR-K-prediction/Data/Aggregated_Data/';
    
    
%     rawBaseDir = 'I:\My Drive\USGS Project\USGS Data\';
%     baseDir = 'I:\My Drive\USGS Project\NMR-K-prediction/Data/Aggregated_Data/';
    
    if strcmp(baseName,'Site1-WellG5')
        site = 'Site1-WellG5';
        name = 'G5_W1_tr5_20x_16p5_up_F1n2_wRIN_wRFI_Reg50_Va1';
        nmrName = name;
    elseif  strcmp(baseName,'Site1-WellG5above')
        site = 'Site1-WellG5';
        name = 'G5_W1_tr5_20x_16p5_up_F1n2_wRIN_wRFI_Reg50_Va1';
        nmrName = 'G5_W1_tr5_20x_16p5_up_F1n2_wRIN_wRFI_Reg50_Va1_above';
    elseif  strcmp(baseName,'Site1-WellG5below')
        site = 'Site1-WellG5';
        name = 'G5_W1_tr5_20x_16p5_up_F1n2_wRIN_wRFI_Reg50_Va1';
        nmrName = 'G5_W1_tr5_20x_16p5_up_F1n2_wRIN_wRFI_Reg50_Va1_below';
    elseif strcmp(baseName,'Site1-WellG6')
        site = 'Site1-WellG6'
        name = 'G6_W2_tr5_20x_16p75_up_F_wRIN_wRFI_reg50_Va1';
        nmrName = name;
    elseif strcmp(baseName,'Site1-WellG6above')
        site = 'Site1-WellG6';
        name = 'G6_W2_tr5_20x_16p75_up_F_wRIN_wRFI_reg50_Va1';
        nmrName = 'G6_W2_tr5_20x_16p75_up_F_wRIN_wRFI_reg50_Va1_above'
    elseif strcmp(baseName,'Site1-WellG6below')
        site = 'Site1-WellG6';
        name = 'G6_W2_ts5_20x_16p75_up_F_wRIN_wRFI_reg50_Va1';
        nmrName = 'G6_W2_tr5_20x_16p75_up_F_wRIN_wRFI_reg50_Va1_below'
    elseif strcmp(baseName,'Site2-WellPN1')
        site = 'Site2-WellPN1';
        name = 'Pl_W1_Tr5_20x_MPp75aLS_F1n2_wRIN_wRFI_Reg50_Va1';
        nmrName = name;
    elseif strcmp(baseName,'Site2-WellPN2')
        site = 'Site2-WellPN2';
        name = 'W2_Tr5_20x_MPp75aLS_Reg50_wRIN_wRFI_Va1';
        nmrName = name;
    else
        nmrName = baseName;
    end

    [d, K, T2ML, phi, z, SumEch, logK, logT2ML, logPhi, SumEch_3s, SumEch_twm, ...
    SumEch_twm_3s] = loadnmrdata2(nmrName);

    in1 = [rawBaseDir site '/' name '/' name '_T2_dist_uniform' '.txt']; 
    in2 = [rawBaseDir site '/' name '/' name '_T2_bins_log10s' '.txt']; 
    in3 = [rawBaseDir site '/' name '/' name '_1Dvectors_uniform' '.txt'];

    T2dist = load(in1);
    T2logbins = load(in2);
    dataVectors = dlmread(in3,'',1,0);
    
    phiVector = dataVectors(:,2);
    
    T2dist(:,1) = T2dist(:,1) - depthOffsets(2); 
       
    T2depths = T2dist(:,1);
    T2data = T2dist(:,2:end);
    
    T2axis = 10.^T2logbins;
    T2axis = T2axis * 10^3;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Find T2 distributions where we have DPP data
    
    for jj = 1:length(z)
       [minVal, index] = min(abs(T2depths - z(jj)));
       filtIndices(jj) = index; 
       filtT2dist(jj,:) = T2data(index,:);        
       filtT2depths(jj) = T2depths(index);
       filtPhi(jj) = phiVector(index);
    end
        
    [sortedK, sortInd] = sort(K);
    sortedT2dist = filtT2dist(sortInd,:);
    sortedPhi = filtPhi(sortInd)';
    T2ML = T2ML(sortInd)';
    z = z(sortInd)';
    
    % Now estimate porosity weighted T2 
    T2linbins = 10.^(T2logbins);

    for jj = 1:length(z)
        sliceT2 = sortedT2dist(jj,:);
        T2W(jj,1) = sum(T2linbins.*sliceT2)./filtPhi(jj);
    end    
    
    T2estimator = T2W;
    
    %%
    for jj = 1:length(T2estimator)
        sliceT2 = sortedT2dist(jj,:);
    % Identify bimodal distributions
       [pks{jj}, locs{jj},widths{jj},proms{jj}] = findpeaks(sliceT2);

       if length(pks{jj}) > 1
          T2peakAvg(jj) =  mean(T2linbins(locs{jj}));
          currentProm = proms{jj};
          currentWidth = widths{jj};
          currentLocs = locs{jj};

          promRatio(jj) = currentProm(1)/currentProm(2);
          widthRatio(jj) = currentWidth(1)/currentWidth(2);

          T2peakMax(jj) = T2linbins(currentLocs(2));
          T2peakMin(jj) = T2linbins(currentLocs(1));
       else
          T2peakAvg(jj) = T2linbins(locs{jj});
          promRatio(jj) = NaN;
          widthRatio(jj) = NaN;
          T2peakMax(jj) = T2linbins(locs{jj});
          T2peakMin(jj) = NaN;
       end
    end
    
    %%
    
     % Try to figure out a range of pore sizes based on bimodal processing
   
    for jj = 1:length(T2estimator)
        if promRatio(jj) > 0.1
            T2_lowerLim(jj) = T2peakMin(jj);
            T2_upperLim(jj) = T2W(jj);
        else
            T2_lowerLim(jj) = T2W(jj);
            T2_upperLim(jj) = T2W(jj);
        end
    end
    
    %%
     SDR_K = @(b,m,n,phi,T2estimator) b.*(phi.^m).*(T2estimator).^n;
     
     b = 0.0035;
     SDR_K_min = SDR_K(b,0,2,phi,T2_lowerLim');
     SDR_K_max = SDR_K(b,0,2,phi,T2_upperLim');
     
     k_est = [SDR_K_min SDR_K_max];
     k_avg = [SDR_K_min SDR_K_max];
     
     k_names = [{'Min T2'} {'Max T2'}];
     
     plotKestKdpp_wdepth(K,k_est,k_avg,z,k_names,2)
     
     

     
end
    