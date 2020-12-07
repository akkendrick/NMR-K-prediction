% Plot basic K, T2ml, water table, phi profiles

clear
%close all

%load bProfileData.mat
%load bestKModels.mat
sites = {'Site1-WellG5','Site1-WellG6','Site2-WellPN1','Site2-WellPN2'};
 
 
gammaNames = {'Site1-G5-Well1-gamma-EMI-bLS.csv','Site1-G6-Well2-gamma-EMI-bLS.csv','Site2-Well1-gamma-EMI-bLS.csv','Site2-Well2-gamma-EMI-bLS.csv'};

%waterTable = [2.0469,2.1248,5.0285,4.7476]; % rel ground surface %NOTE G5/G6 cased below clay,
%   need to use nearby water level from above the clay, for G5 + G6 using
%   water level from well well G2 cased above the New Rome Clay (rel to gs)
waterTable = [2.0469,2.1248,5.0285,4.7476];


 for jj = 1:length(sites)
    baseName = sites{jj}
%      baseDir = '/Volumes/GoogleDrive/My Drive/Stanford/USGS Project/Field Data/USGS Data/';
%      gammaBaseDir = '/Volumes/GoogleDrive/My Drive/Stanford/USGS Project/Field Data/USGS Data/WI_gamma-EMI-bLS_csvfiles';

    baseDir = 'I:\My Drive\Stanford\USGS Project\Field Data\USGS Data\';
    gammaBaseDir = 'I:\My Drive\Stanford\USGS Project\Field Data\USGS Data\WI_gamma-EMI-bLS_csvfiles';
    
    [T2dist,T2logbins,SEdecayTime,SEdecayUniform,SEdecay,oneDVectors,...
        oneDVectorsUniform, nmrName] = loadAllRawNMRdata(baseName);

    in3 = [baseDir baseName '/' strcat(baseName,'_DPP_filt.txt')];
    in7 = [gammaBaseDir '/' gammaNames{jj}];
    
    DPPdat = load(in3); 
    gammaEMIdata = csvread(in7,3,0);
    
%     T2depths = T2dist(:,1);
%     T2dist = T2dist(:,2:end);

    [d, K, T2ML, phi, z, SumEch, logK, logT2ML, logPhi, SumEch_3s, SumEch_twm, ...
    SumEch_twm_3s] = loadnmrdata2(baseName); 

    NMRphi = oneDVectors(:,2);

    Dk = DPPdat(:,2)*1.16e-5; % converts K from m/day to m/s
    z_dk = DPPdat(:,1);

    origT2dist = T2dist;
    
    %%
    zNMR = z';
   
    %%
    
    gammaEMIoffset = [0,0,0.75, 0.75];
    
    if jj <= d
        
%         plotWellProfiles_wLith(K,NMRphi,waterTable(kk),z,T2ML,T2dist,T2logbins,[],{'K DPP'},{'*'},bProfile{kk},...
%             lithStart{kk}, lithEnd{kk},lithUnits{kk},dlubacModel{kk},bestSDRModel{kk},meanT2B)
        plotWellProfiles_wNMRprofile_wGamma_EMI(K,NMRphi,waterTable(jj),z,T2ML,T2dist,T2logbins,gammaEMIdata,gammaEMIoffset(jj), DPPdat)
           

    else
%         plotWellProfiles(K,NMRphi,waterTable(kk),z,T2ML,T2dist,T2logbins,[],{'K DPP'},{'*'},bProfile{kk},...
%            dlubacModel{kk},bestSDRModel{kk},meanT2B)
        plotWellProfiles_wNMRprofile_wGamma_EMI(K,NMRphi,waterTable(jj),z,T2ML,T2dist,T2logbins,gammaEMIdata,gammaEMIoffset(jj),DPPdat)

    end
    
%     fileName = strcat(baseName,'_profile.fig');
%     savefig(fileName)
%     
%     fileName = strcat(baseName,'_profile.png');
%     print('-dpng','-r300',fileName)
%     
%     fileName = strcat(baseName,'_profile.svg');
%     print('-dsvg','-r300',fileName)
    
    %title(baseName)
    
 end
    