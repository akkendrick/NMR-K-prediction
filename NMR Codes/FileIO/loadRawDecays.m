function [decayCurves,decayTime] = loadRawDecays(siteName)
%LOADRAWDECAYS 
% load raw NMR data for processing into decay curves

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% replace with path to data
rawBaseDir = 'E:\Dropbox\Research\Lab Files\Field Data\USGS Data\';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

wisc_sites = {'Site1-WellG5','Site1-WellG6','Site2-WellPN1','Site2-WellPN2'};

if sum(strcmp(siteName,wisc_sites) == 1)
    if strcmp(siteName,'Site1-WellG5')
        site = 'Site1-WellG5';
        name = 'G5_W1_tr5_20x_16p5_up_F1n2_wRIN_wRFI_Reg50_Va1';
        nmrName = site;
        offset = 0.95;
        %offset = 3.12422224
    elseif strcmp(siteName,'Site1-WellG6')
        site = 'Site1-WellG6';
        name = 'G6_W2_tr5_20x_16p75_up_F_wRIN_wRFI_reg50_Va1';
        nmrName = site;
        offset = 0.75;
    elseif strcmp(siteName,'Site2-WellPN1')
        site = 'Site2-WellPN1';
        name = 'Pl_W1_Tr5_20x_MPp75aLS_F1n2_wRIN_wRFI_Reg50_Va1';
        nmrName = site;
        offset = 0.75;
    elseif strcmp(siteName,'Site2-WellPN2')
        site = 'Site2-WellPN2';
        name = 'W2_Tr5_20x_MPp75aLS_Reg50_wRIN_wRFI_Va1';
        nmrName = site;
        offset = 0.75;
    else
        disp('Not a good filename')
        return
    end
        in1 = [rawBaseDir '/' site '/' site '_SE_decay' '.txt'];
        in2 = [rawBaseDir '/' site '/' site '_SE_decay_time' '.txt'];

        decayCurves = load(in1);
        decayTime = load(in2);
        
else
    disp('Not a good filename')
    return    

end

end

