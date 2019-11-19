function [gammaEMIdepth, gamma, EMI] = loadGammaEMIData(siteName)

    kanWash_sites = {'dpnmr_larned_east','dpnmr_larned_lwph','dpnmr_larned_west',...
    'dpnmr_leque_east','dpnmr_leque_west','dpnmrA11','dpnmrA12',...
    'dpnmrC1S','dpnmrC1SE','dpnmrC1SW'};

    wisc_sites = {'Site1-WellG5','Site1-WellG6','Site2-WellPN1','Site2-WellPN2'};

    if sum(strcmp(siteName,kanWash_sites) == 1)
        disp('No gamma/EMI data')
        return
    elseif sum(strcmp(siteName,wisc_sites) == 1)
        gammaBaseDir = 'I:\My Drive\USGS Project\USGS Data\WI_gamma-EMI-bLS_csvfiles';
        %gammaBaseDir = '/Volumes/GoogleDrive/My Drive/USGS Project/USGS Data/WI_gamma-EMI-bLS_csvfiles';

        
        
        if strcmp(siteName,'Site1-WellG5')
            site = 'Site1-WellG5';
            name = 'G5_W1_tr5_20x_16p5_up_F1n2_wRIN_wRFI_Reg50_Va1';
            gammaName = 'Site1-G5-Well1-gamma-EMI-bLS.csv';
            nmrName = name;
            offset = 0.95;
        elseif strcmp(siteName,'Site1-WellG6')
            site = 'Site1-WellG6';
            name = 'G6_W2_tr5_20x_16p75_up_F_wRIN_wRFI_reg50_Va1';
            gammaName = 'Site1-G6-Well2-gamma-EMI-bLS.csv';
            nmrName = name;
            offset = 0.75;
        elseif strcmp(siteName,'Site2-WellPN1')
            site = 'Site2-WellPN1';
            name = 'Pl_W1_Tr5_20x_MPp75aLS_F1n2_wRIN_wRFI_Reg50_Va1';
            gammaName = 'Site2-Well1-gamma-EMI-bLS.csv';
            nmrName = name;
            offset = 0.75;
        elseif strcmp(siteName,'Site2-WellPN2')
            site = 'Site2-WellPN2';
            name = 'W2_Tr5_20x_MPp75aLS_Reg50_wRIN_wRFI_Va1';
            gammaName = 'Site2-Well2-gamma-EMI-bLS.csv';
            nmrName = name;
            offset = 0.75;
        else
            disp('Not a good filename')
            return
        end
        
        in7 = [gammaBaseDir '/' gammaName]
        gammaEMIdata = csvread(in7,3,0);

        gammaEMIdepth = gammaEMIdata(:,1);
        gamma = gammaEMIdata(:,2);
        EMI = gammaEMIdata(:,3);
        
    end


end

