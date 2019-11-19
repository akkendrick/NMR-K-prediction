% Compare T2 distributions at similar depths

sites = {'Site1-WellG5','Site1-WellG6'}

%waterTable = [2.0469,2.1248,5.0285,4.7476]; % rel ground surface %NOTE G5/G6 cased below clay,
%   need to use nearby water level from above the clay, for G5 + G6 using
%   water level from well well G2 cased above the New Rome Clay (rel to gs)
waterTable = [1.181,1.181,5.0285,4.7476];


%waterTable = [3.004,2.963,5.727,5.408]; % rel top of casing
depthOffsets = [0.75,0.95,0.75,0.75];

for kk = 1:length(sites)
     baseName = sites{kk}
     baseDir = '/Volumes/GoogleDrive/My Drive/USGS Project/USGS Data/';
     gammaBaseDir = '/Volumes/GoogleDrive/My Drive/USGS Project/USGS Data/WI_gamma-EMI-bLS_csvfiles';

%     baseDir = 'I:\My Drive\USGS Project\USGS Data\';
%     gammaBaseDir = 'I:\My Drive\USGS Project\USGS Data\WI_gamma-EMI-bLS_csvfiles';

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
        site = 'Site1-WellG6';
        name = 'G6_W2_tr5_20x_16p75_up_F_wRIN_wRFI_reg50_Va1';
        nmrName = name;
    elseif strcmp(baseName,'Site1-WellG6above')
        site = 'Site1-WellG6';
        name = 'G6_W2_tr5_20x_16p75_up_F_wRIN_wRFI_reg50_Va1';
        nmrName = 'G6_W2_tr5_20x_16p75_up_F_wRIN_wRFI_reg50_Va1_above'
    elseif strcmp(baseName,'Site1-WellG6below')
        site = 'Site1-WellG6';
        name = 'G6_W2_tr5_20x_16p75_up_F_wRIN_wRFI_reg50_Va1';
        nmrName = 'G6_W2_tr5_20x_16p75_up_F_wRIN_wRFI_reg50_Va1_below'
    elseif strcmp(baseName,'Site2-WellPN1')
        site = 'Site2-WellPN1';
        name = 'Pl_W1_Tr5_20x_MPp75aLS_F1n2_wRIN_wRFI_Reg50_Va1';
        nmrName = name;
    elseif strcmp(baseName,'Site2-WellPN2')
        site = 'Site2-WellPN2';
        name = 'W2_Tr5_20x_MPp75aLS_Reg50_wRIN_wRFI_Va1';
        nmrName = name;
    end

    in1 = [baseDir site '/' name '/' name '_SE_decay' '.txt']; 
    in2 = [baseDir site '/' name '/' name '_1Dvectors.txt'];
    in3 = [baseDir site '/' name '/' strcat(site,'_DPP.txt')];
    in4 = [baseDir site '/' name '/' name '_SE_decay_time.txt'];
    in5 = [baseDir site '/' name '/' name '_T2_dist.txt'];
    in6 = [baseDir site '/' name '/' name '_T2_bins_log10s.txt'];
    
    decaycurve = load(in1); 
    dparam = dlmread(in2,'\t',1,0); 
    DPPdat = load(in3); 
    t = load(in4);
    
    T2dist = load(in5);
    T2dist(:,1) = T2dist(:,1) - depthOffsets(kk);
    
    T2dists{kk} = T2dist;
    
    T2depths = T2dist(:,1);
    T2logbins = load(in6);
    
    % set needed variables
    NMRphi = dparam(:,2);
    S = decaycurve(:,2:end); 

    %Dk = olddat(:,4)*1.16e-5; 
    %z_dk = olddat(:,1); 

    [d, K, T2ML, phi, z, SumEch, logK, logT2ML, logPhi, SumEch_3s, SumEch_twm, ...
    SumEch_twm_3s] = loadnmrdata2(nmrName); 

    Dk = DPPdat(:,2)*1.16e-5; % converts K from m/day to m/s
    z_dk = DPPdat(:,1);

    origT2dist = T2dist;
    
    %%
    zNMR{kk} = T2depths;
    
    T2linbins = 10.^T2logbins;
    
end

% Organize data for comparision
zNMR_filt{1} = zNMR{1}(39:63);
zNMR_filt{2} = zNMR{2}(40:64);

T2dists_filt{1} = T2dists{1}(39:63,:);
T2dists_filt{2} = T2dists{2}(40:64,:);

for jj = 1:24
    for kk = 1:length(sites)

        currentT2dist = T2dists_filt{kk}(:,2:end);
        currentZ = zNMR_filt{kk};
        
        figure(jj)
        hold on
        
        plot(T2linbins, currentT2dist(jj,:), 'LineWidth',2)

        grid on
        box on
      
        set(gca,'XScale','log')
        
        legendString{kk} = strcat('z= ', string(currentZ(jj)),' ',sites{kk});
    
    end
    
    legend(legendString{:})
    
end