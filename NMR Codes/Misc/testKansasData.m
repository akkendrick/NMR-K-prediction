%baseDir = 'I:\My Drive\USGS Project\Kansas_Wash_Data\';
baseDir = '/Volumes/GoogleDrive/My Drive/USGS Project/Kansas_Wash_Data/';

name = 'dpnmr_larned_east.txt';
site = 'dpnmr_larned_east';
nmrName = 'dpnmr_larned_east';

%in3 = [baseDir site '/' strcat(site,'_DPP_filt.txt')];
in5 = [baseDir site '/' strcat(site, '_T2_dist.txt')];
%in6 = [baseDir site '/' name '/' name '_T2_bins_log10s.txt'];

%DPPdat = load(in3); 
T2dist = load(in5);
%T2logbins = load(in6);

[d, K, T2ML, phi, z, SumEch, logK, logT2ML, logPhi, SumEch_3s, SumEch_twm, ...
SumEch_twm_3s] = loadnmrdata2(nmrName); 