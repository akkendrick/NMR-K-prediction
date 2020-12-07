clear
close all
%name = 'all_data'
%name = 'A1'
%name = 'C1'
%name = 'dpnmr_larned_east'
%name = 'dpnmr_larned_west'
%name = 'dpnmr_larned_lwph'
%name = 'dpnmr_leque_east'
%name = 'dpnmr_leque_all'

%name = 'dpnmr_leque_west'

%name = 'wisc_all';

name = 'Site1-WellG5';

%name = 'dpnmr_leque_west'

%baseDir = '/Volumes/GoogleDrive/My Drive/Stanford/USGS Project/Field Data/USGS Data/';
baseDir = 'I:\My Drive\Stanford\USGS Project\Field Data\USGS Data\';

%baseDir = '/Volumes/GoogleDrive/My Drive/Stanford/USGS Project/Kansas_Wash_Data/';

site = name;
if strcmp(site,'Site1-WellG5')
    name = 'G5_W1_tr5_20x_16p5_up_F1n2_wRIN_wRFI_Reg50_Va1';
    nmrName = name;
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
elseif strcmp(site,'Site1-WellG6above')
    site = 'Site1-WellG6';
    name = 'G6_W2_tr5_20x_16p75_up_F_wRIN_wRFI_reg50_Va1';
    nmrName = 'G6_W2_tr5_20x_16p75_up_F_wRIN_wRFI_reg50_Va1_above'
elseif strcmp(site,'Site1-WellG6below')
    site = 'Site1-WellG6';
    name = 'G6_W2_tr5_20x_16p75_up_F_wRIN_wRFI_reg50_Va1';
    nmrName = 'G6_W2_tr5_20x_16p75_up_F_wRIN_wRFI_reg50_Va1_below'
elseif strcmp(site,'Site2-WellPN1')
    name = 'Pl_W1_Tr5_20x_MPp75aLS_F1n2_wRIN_wRFI_Reg50_Va1';
    nmrName = name;
elseif strcmp(site,'Site2-WellPN2')
    name = 'W2_Tr5_20x_MPp75aLS_Reg50_wRIN_wRFI_Va1';
    nmrName = name;
else
    nmrName = name;
end

in1 = [baseDir site '/' name '/' name '_T2_dist' '.txt']; 
in2 = [baseDir site '/' name '/' name '_T2_bins_log10s' '.txt']; 

%in1 = [baseDir site '/' name '_T2_dist' '.txt']; 
%in2 = [baseDir site '/' name '_T2_bins_log10s' '.txt']; 



T2dist = load(in1); 
T2logbins = load(in2);


[d, K, T2ML, phi, z, SumEch, logK, logT2ML, logPhi, SumEch_3s, SumEch_twm, ...
    SumEch_twm_3s] = loadnmrdata2(nmrName); 

m = 1;
[KGM_lk, bestTau, bestRho, r] = grid_search_kgm(logT2ML, phi, logK, m);

KGM_k = 10.^KGM_lk;

display('RMSE Error')
bestError = computeError(K, KGM_k)

plotKwithDepth(K,z,T2dist,T2logbins,KGM_k,{'DPP','KGM'},{'+'})

title(name)

%[bestTau, bestRho, r] = grid_search_KGM_akk(T2ML, phi, K);