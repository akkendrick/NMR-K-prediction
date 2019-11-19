% Compute mean pore size from K and NMR logging data

sites = {'wisc_all'}
legendNames = {'Wisc'}

for kk = 1:length(sites)
    %close all
    baseName = sites{kk}
    baseDir = '/Volumes/GoogleDrive/My Drive/USGS Project/USGS Data/';
    baseDir = '/Volumes/GoogleDrive/My Drive/USGS Project/NMR-K-prediction/Data/Aggregated_Data/';

    %baseDir = 'I:\My Drive\USGS Project\USGS Data\';
    
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

    SDR_K = @(b,m,n,phi,T2ML) (b.*(phi.^m).*(T2ML).^n);
    SDR_b = @(K,m,n,phi,T2ML) K./((phi.^m).*(T2ML).^n);
    
    SOE_K = @(C,SOE,n) C.*(SOE.^n);
    SOE_C = @(K,SOE,n) K./(SOE.^n);
    
    m = 0;
    n = 2;
    
    bProfile{kk} = SDR_b(K,m,n,phi,T2ML);
    cProfile{kk} = SOE_C(K, SumEch, n);
    
    %%
    % Sort data and correct large T2ML
    
    [sortedLogT2ML, sortInd] = sort(logT2ML);
    sortedLogK = logK(sortInd,:);
    
    %%
    % ESTIMATE MEAN PORE SIZE
    Kcm = K *10^2; % m/s -> cm/s
    
    voidRatio = phi ./ (1-phi);
    refVoidRatio = 1;
    
    beta = 4;
    refK = Kcm.*(voidRatio./refVoidRatio).^(-beta);
    
    logSpecSurf = (log10(refK)+5)./(-2);
    specSurf = 10.^(logSpecSurf); % in m^2/g
           
    rho_m = 2650; % kg/m^3 for quartz, seems representative
    
    meanPoreSize = 2*(voidRatio)./(specSurf.*rho_m); %in m
    meanPoreSizeMM = meanPoreSize.*10^3;
    
    figure(1)
    grid on
    hold on
    
    scatter(meanPoreSize,Kcm,40,'Filled')
    
    set(gca,'XScale','log')
    set(gca,'YScale','log')
    
    xlabel('Estimated mean pore size (mm)')
    ylabel('Hydraulic Conductivity (cm/s)')
    
    box on
    
    set(gca,'FontSize',12)
    
    figure(2)
    grid on
    hold on
    
    scatter(meanPoreSize, T2ML.^2,40,'Filled')
    set(gca,'XScale','log')
    set(gca,'YScale','log')
    
    xlabel('Estimated mean pore size (mm)')
    ylabel('T_{2ML}^2 (s)')
    
    box on
    
    set(gca,'FontSize',12)
    
      
    % Attempt to estimate pore radius from Godefroy
    D = @(temp) (1.0413+0.039828*temp+0.00040318*temp^2)*10^-9;
    T2B = 1;
    alpha = 3;
    %bParam = (T2ML.*T2B.*alpha)./(T2B - T2ML);
    
    SeeversT2 = (T2ML.^(-1) - T2B^(-1)).^(-1);
    
    %bParam = (T2ML.*T2B.*alpha)./(T2B - T2ML);
    bParam = (SeeversT2.*T2B.*alpha)./(T2B-SeeversT2);
    
    
    
    rho = 44*10^-6; %m/s
    rhoMod = 2;
    rho = rho * rhoMod;
    
    Dtemp = D(15);
    
    poreSize = (sqrt(Dtemp.*(2.*bParam.*rho^2+Dtemp)) - Dtemp)./rho;
    
    figure(2)
    scatter(poreSize*10^3, T2ML.^2,30,'Filled')

    figure(1)
    scatter(poreSize*10^3, Kcm, 30, 'Filled') 
    
    % Attempt to estimate an r from rho S/V
    rEst = (3*T2ML*rho);
    dEst = 2*rEst;
    
    figure(2)
    scatter(dEst*10^3, T2ML.^2,30,'Filled')
    
    
    
end
    
    
