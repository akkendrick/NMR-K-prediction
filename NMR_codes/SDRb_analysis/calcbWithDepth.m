% Directly estimate b as a function of depth

close all
clear

siteList = [{'Site1-WellG6'} {'Site1-WellG6above'} {'Site1-WellG6below'}...
    {'Site1-WellG5'} {'Site1-WellG5above'} {'Site1-WellG5below'} {'Site2-WellPN1'} {'Site2-WellPN2'}];

baseDir = '/Volumes/GoogleDrive/My Drive/USGS Project/USGS Data/';
%baseDir = 'I:\My Drive\USGS Project\USGS Data\';

load SDR_bestFit_table_m0_n2.mat

for kk = 1:length(siteList)
    
    site = siteList{kk};

    if strcmp(site,'Site1-WellG5')
        name = 'G5_W1_tr5_20x_16p5_up_F1n2_wRIN_wRFI_Reg50_Va1';
        nmrName = name;
        saveName = 'G5';

    elseif  strcmp(site,'Site1-WellG5above')
        site = 'Site1-WellG5';
        name = 'G5_W1_tr5_20x_16p5_up_F1n2_wRIN_wRFI_Reg50_Va1';
        nmrName = 'G5_W1_tr5_20x_16p5_up_F1n2_wRIN_wRFI_Reg50_Va1_above';

    elseif  strcmp(site,'Site1-WellG5below')
        site = 'Site1-WellG5';
        name = 'G5_W1_tr5_20x_16p5_up_F1n2_wRIN_wRFI_Reg50_Va1';
        nmrName = 'G5_W1_tr5_20x_16p5_up_F1n2_wRIN_wRFI_Reg50_Va1_below';

    elseif strcmp(site,'Site1-WellG6')
        name = 'G6_W2_tr5_20x_16p75_up_F_wRIN_wRFI_reg50_Va1';
        nmrName = name;
        saveName = 'G6';

    elseif strcmp(site,'Site1-WellG6above')
        site = 'Site1-WellG6';
        name = 'G6_W2_tr5_20x_16p75_up_F_wRIN_wRFI_reg50_Va1';
        nmrName = 'G6_W2_tr5_20x_16p75_up_F_wRIN_wRFI_reg50_Va1_above';

    elseif strcmp(site,'Site1-WellG6below')
        site = 'Site1-WellG6';
        name = 'G6_W2_tr5_20x_16p75_up_F_wRIN_wRFI_reg50_Va1';
        nmrName = 'G6_W2_tr5_20x_16p75_up_F_wRIN_wRFI_reg50_Va1_below';

    elseif strcmp(site,'Site2-WellPN1')
        name = 'Pl_W1_Tr5_20x_MPp75aLS_F1n2_wRIN_wRFI_Reg50_Va1';
        nmrName = name;
        saveName = 'PN1';

    elseif strcmp(site,'Site2-WellPN2')
        name = 'W2_Tr5_20x_MPp75aLS_Reg50_wRIN_wRFI_Va1';
        nmrName = name;
        saveName = 'PN2';

    else
        name = site;
        nmrName = site;
        saveName = 'wisc_all';
    end
    
    %if figureson == 1
    in1 = [baseDir site '/' name '/' name '_T2_dist' '.txt']; 
    in2 = [baseDir site '/' name '/' name '_T2_bins_log10s' '.txt']; 

    T2dist = load(in1); 
    T2logbins = load(in2);
    
    SDR_K = @(b,m,n,phi,T2ML) b.*(phi.^m).*(T2ML).^n;
    SDR_b = @(K,m,n,phi,T2ML) K./((phi.^m).*(T2ML).^n);

    [d, K, T2ML, phi, z, SumEch, logK, logT2ML, logPhi, SumEch_3s, SumEch_twm, ...
    SumEch_twm_3s] = loadnmrdata2(nmrName); 

    logSumEch = log10(SumEch);

    bestFitMatrix = zeros(3,3);
    totalErrorEstimate = zeros(1,3);
    
    current_n = SDR_Bootstrap(2,kk);
    current_m = SDR_Bootstrap(3,kk);
    
    %currrent_m = 0;
    %current_n = 2;
    
    bProfileTemp{:,kk} = SDR_b(K,current_m,current_n,phi,T2ML);
    
    kProfileTemp{:,kk} = SDR_K(bProfileTemp{:,kk},current_m,current_n,phi,T2ML);
     
end

bProfile{:,1} = vertcat(bProfileTemp{:,6},bProfileTemp{:,5});
bProfile{:,2} = vertcat(bProfileTemp{:,3},bProfileTemp{:,2});
bProfile{:,3} = bProfileTemp{:,7};
bProfile{:,4} = bProfileTemp{:,8};

kProfile{:,1} = vertcat(kProfileTemp{:,6},kProfileTemp{:,5});
kProfile{:,2} = vertcat(kProfileTemp{:,3},kProfileTemp{:,2});
kProfile{:,3} = kProfileTemp{:,7};
kProfile{:,4} = kProfileTemp{:,8};

save('bProfileData.mat','bProfile','siteList','kProfile')


