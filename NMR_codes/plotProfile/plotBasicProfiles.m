% Plot basic K, T2ml, water table, phi profiles

clear
%close all

load bProfileData.mat
load bestKModels.mat
load('sonicCoreT2B.mat','sonicCoreT2BData')

sites = {'Site1-WellG5','Site1-WellG6',...
   'Site2-WellPN1','Site2-WellPN2'};
 
% sites = {'Site2-WellPN1'}
 
gammaNames = {'Site1-G5-Well1-gamma-EMI-bLS.csv','Site1-G6-Well2-gamma-EMI-bLS.csv',...
    'Site2-Well1-gamma-EMI-bLS.csv','Site2-Well2-gamma-EMI-bLS.csv'};

%waterTable = [2.0469,2.1248,5.0285,4.7476]; % rel ground surface %NOTE G5/G6 cased below clay,
%   need to use nearby water level from above the clay, for G5 + G6 using
%   water level from well well G2 cased above the New Rome Clay (rel to gs)
waterTable = [1.181,1.181,5.0285,4.7476];


%waterTable = [3.004,2.963,5.727,5.408]; % rel top of casing
depthOffsets = [0.75,0.95,0.75,0.75];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filter sonic core logging data (ignore mysterious points with short relax
% times likely due to oxidized Fe2+ after sample collection

filteredT2Bdata = sonicCoreT2BData.T2Bpeak(sonicCoreT2BData.T2Bpeak > 1000);
meanT2B = mean(filteredT2Bdata);
meanT2B = meanT2B * 10^-3; % Put data in seconds

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute lithology logs from NMR data
% clay = 1 sand = 2 gravel = 3
% G5start = [0,9.95,11.2,11.7,14.7];
% G5end = [9.95,11.2,11.7,14.7,16.7];
% G5lith = [2,1,2,1,2];
% 
% G6start = [0, 12, 14.25];
% G6end = [12,14.25,16.7];
% G6lith = [2,1,2];
% 
% PN1start = [0];
% PN1end = [15.82];
% PN1lith = [2];
% 
% PN2start = [0];
% PN2end = [17.7];
% PN2lith = [2];
% 
% lithStart = {G5start, G6start, PN1start, PN2start};
% lithEnd = {G5end, G6end, PN1end, PN2end};
% lithUnits = {G5lith, G6lith, PN1lith, PN2lith};

% Using Dave's lithology
G5start = [0,35,47]./3.281;
G5end = [35,47,65]./3.281;
G5lith = [2,1,2];

G6start = [0,35,47]./3.281;
G6end = [35,47,65]./3.281;
G6lith = [2,1,2];

lithStart = {G5start,G6start};
lithEnd = {G5end,G6end};
lithUnits = {G5lith, G6lith};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    in7 = [gammaBaseDir '/' gammaNames{kk}];
    
    gammaEMIdata = csvread(in7,3,0);

    load('sonicCoreT2B.mat','sonicCoreT2BData')
    T2B_depth = sonicCoreT2BData.Depthm;
    T2B_depth = flipud(T2B_depth);  
    
    T2B_peak = sonicCoreT2BData.T2Bpeak*(10^-3);
    T2B_peak = flipud(T2B_peak); 
    
    decaycurve = load(in1); 
    dparam = dlmread(in2,'\t',1,0); 
    DPPdat = load(in3); 
    t = load(in4);
    
    T2dist = load(in5);
    T2dist(:,1) = T2dist(:,1) - depthOffsets(kk);
    
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
    zNMR = z';
    
    
        
%     %%
%     % Quickly plot a T2 slice
%     
%     % Pick out representative slices in G5
%     % Short K slice (z = 3.565 use 3.45)
%     % Large K slice (z = 7.55 use 7.45)
%     % Clay slice (z = 13.45)
%     
%     shortKSlice = T2dist(54,2:end);
%     largeKSlice = T2dist(38,2:end);
%     clayKSlice = T2dist(14,2:end);
%     
%     sliceT2 = T2dist(40,2:end);
%     T2linbins = 10.^T2logbins;
%     
%     T2ML_shortSlice = 10.^(sum(T2logbins.*shortKSlice)./dparam(54,2));
%     [minVal, ind] = min(abs(T2linbins - T2ML_shortSlice));
%     
%     T2ML_shortAmp = shortKSlice(ind);
%     
%     T2ML_largeSlice = 10.^(sum(T2logbins.*largeKSlice)./dparam(38,2));
%     [minVal, ind] = min(abs(T2linbins - T2ML_largeSlice));
%     
%     T2ML_largeAmp = largeKSlice(ind);
%     
%     T2ML_claySlice = 10.^(sum(T2logbins.*clayKSlice)./dparam(14,2));
%     [minVal, ind] = min(abs(T2linbins - T2ML_claySlice));
%     
%     T2ML_clayAmp = clayKSlice(ind);
%     
%     
%     figure(10)
%     hold on
%     plot(T2linbins, largeKSlice,'LineWidth',3)
%     plot(T2linbins, shortKSlice,'LineWidth',3)
%     plot(T2linbins, clayKSlice,'LineWidth',3)
% 
%     plot(T2ML_largeSlice, T2ML_largeAmp,'r*','MarkerSize',15)
%     plot(T2ML_shortSlice, T2ML_shortAmp,'r*','MarkerSize',15)
%     plot(T2ML_claySlice, T2ML_clayAmp,'r*','MarkerSize',15)
% 
%     
%     grid on
%     box on
%     set(gca,'XScale','log')
%     
%     xlabel('Relaxation Time T_2 (s)')
%     ylabel('Amplitude')
%     xlim([0.0007875, 5])
%     
%     legend('High K Slice','Low K Slice','Clay Slice','T_{2ML}')
%     
%     set(gca,'FontSize',14)
%     
    
    %%
    
    if kk <= 2
        
%         plotWellProfiles_wLith(K,NMRphi,waterTable(kk),z,T2ML,T2dist,T2logbins,[],{'K DPP'},{'*'},bProfile{kk},...
%             lithStart{kk}, lithEnd{kk},lithUnits{kk},dlubacModel{kk},bestSDRModel{kk},meanT2B)
        plotWellProfiles_wGamma_EMI(K,NMRphi,waterTable(kk),z,T2ML,T2dist,T2logbins,gammaEMIdata, DPPdat)
           

    else
%         plotWellProfiles(K,NMRphi,waterTable(kk),z,T2ML,T2dist,T2logbins,[],{'K DPP'},{'*'},bProfile{kk},...
%            dlubacModel{kk},bestSDRModel{kk},meanT2B)
        plotWellProfiles_wGamma_EMI(K,NMRphi,waterTable(kk),z,T2ML,T2dist,T2logbins,gammaEMIdata,DPPdat)

    end
    
    fileName = strcat(baseName,'_profile.fig');
    savefig(fileName)
    
    fileName = strcat(baseName,'_profile.png');
    print('-dpng','-r300',fileName)
    
    fileName = strcat(baseName,'_profile.svg');
    print('-dsvg','-r300',fileName)
    
    %title(baseName)
    
 end
    