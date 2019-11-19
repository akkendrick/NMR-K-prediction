 % Compare all permeabililty models
clear
close all
 
 sites = {'Site1-WellG5above','Site1-WellG5below','Site1-WellG6above',...
     'Site1-WellG6below','Site2-WellPN1','Site2-WellPN2'}
 
 %sites = {'Site2-WellPN1'}
 
%  SDR_n = [2,0,2,0,1,1];
%  SDR_m = [1,0,1,0,1,1];
%  SDR_b = [0.008,2.75E-04,0.008,2.75E-04,3.8E-03,3.8E-03]; 
%  

 SDR_n = [2,2,2,2,2,2];
 SDR_m = [0,0,0,0,0,0];
 SDR_b = [6.232e-05,6.232e-05,0.0001434,0.0001434,0.0001309,7.824e-05]; 

 cutoffs = [174,800,324,246,244,266]*10^-3;
 TC_n = [1,1,1,1,1,1];
 TC_m = [2,2,2,2,2,2];
 %TC_c = [6.93E-08,3.19E-10,6.93E-08,3.19E-10,0.0122,0.0122];
 TC_c = [0.0036,0.1274,0.0078,0.0104,0.0229,0.0161];
 
%  Seevers_m = [1,1,1,1,1,1];
%  Seevers_n = [2,0,2,0,1,1];
%  Seevers_b = [0.00724,0.000748,0.00724,0.000748,0.00355,0.00355];

 Seevers_m = [1,0,1,0,1,1];
 Seevers_n = [1,1,1,1,1,1];
 Seevers_b = [0.0010,0.0065,0.0010,0.0065,0.0034,0.0034];
 
 KGM_tau = [2.1135,2.1135,1,1.9724,2.6915,3.1989];
 KGM_rho = [6.95E-05,1,1.91E-05,1,1,1];
 KGM_m = [0,0,0,0,0,0];
 
 SOE_n = [1,1,1,1,1,1];
 SOE_b = [0.0066,0.0067,0.0064,0.0063,0.0187,0.0105];
 
 SDR_K = @(b,m,n,phi,T2ML) b.*(phi.^m).*(T2ML).^n;
 TC_K = @(c,m,n,phi,indexQuotient) c.*(phi).^m.*(indexQuotient).^n;
 SOE_K = @(b,n,SOE) b.*(SOE).^n;
 Seevers_K = @(b,m,n,T2ML,T2B,phi) b.*(phi).^m.*(T2ML.^(-1) - T2B.^(-1)).^(-1);
  
 
Temp = 20;  % temperature in degress C 
rho = @(Tt) 1000*(1 - ((Tt+288.94)./(508929*(Tt+68.12))).*(Tt-3.98).^2); % kg/m^3
eta = @(Tt) 0.0013 - 1.7e-5*Tt;         % Pa -s
%Tb = @(Tt) 3.3 + 0.044*(Tt - 35);       % seconds

Tb = @(Tt) 200;
D = @(Tt) (1.0413 + 0.039828*Tt + 0.00040318*Tt.^2).*1e-9;  % m^2/s 
g = 9.8;    %m/s^2
tort = 1/(1.5^2); 
t1 = (rho(Temp)*g)/(8*eta(Temp)); % 

num2 = @(T2) 4*D(Temp)*Tb(Temp)*T2;
denom2 = @(T2) Tb(Temp) - T2; 
    
f12 = @(rho) (D(Temp)./rho);  
SQterm = @(rho,T2) sqrt(f12(rho).^2 + (num2(T2)./denom2(T2))); 

KGM_lK = @(rho,tau,m,lphi,T2) log10(1/tau^2) + log10(t1) + m*lphi + 2*log10(SQterm(rho,T2)-f12(rho)); 

errorEstimates = zeros(length(sites),4);

for kk = 1:length(sites)
    baseName = sites{kk}
    baseDir = '/Volumes/GoogleDrive/My Drive/USGS Project/USGS Data/';

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
    in3 = [baseDir site '/' name '/' strcat(site,'_DPP_filt.txt')];
    in4 = [baseDir site '/' name '/' name '_SE_decay_time.txt'];
    in5 = [baseDir site '/' name '/' name '_T2_dist.txt'];
    in6 = [baseDir site '/' name '/' name '_T2_bins_log10s.txt'];

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
    T2logbins = load(in6);
    % set needed variables
    z = dparam(:,1); 

    S = decaycurve(:,2:end); 

    %Dk = olddat(:,4)*1.16e-5; 
    %z_dk = olddat(:,1); 

    [d, K, T2ML, phi, z, SumEch, logK, logT2ML, logPhi, SumEch_3s, SumEch_twm, ...
    SumEch_twm_3s] = loadnmrdata2(nmrName); 

    Dk = DPPdat(:,2)*1.16e-5; % converts K from m/day to m/s
    z_dk = DPPdat(:,1);

    origT2dist = T2dist;

    zT2dist = T2dist(:,1);
    T2dist = T2dist(:,2:end);
    T2linbins = 10.^T2logbins;

    if strcmp(baseName,'Site1-WellG6above')
        % Take depths shallower than 11m
        z_dk = z_dk(z_dk < 11);
        Dk = Dk(z_dk < 11)

        T2dist = T2dist(zT2dist < 11,:);
        zT2dist = zT2dist(zT2dist < 11,:);
    elseif strcmp(baseName,'Site1-WellG5above')
        % Take depths shallower than 11m
        z_dk = z_dk(z_dk < 11);
        Dk = Dk(z_dk < 11);

        T2dist = T2dist(zT2dist < 11,:);
        zT2dist = zT2dist(zT2dist < 11,:);
    elseif strcmp(baseName,'Site1-WellG6below')
        z_dk = z_dk(z_dk > 11);
        Dk = Dk(z_dk > 11);

        T2dist = T2dist(zT2dist > 11,:);
        zT2dist = zT2dist(zT2dist > 11,:);
    elseif strcmp(baseName,'Site1-WellG5below')
        z_dk = z_dk(z_dk > 11);
        Dk = Dk(z_dk > 11);

        T2dist = T2dist(zT2dist > 11,:);
        zT2dist = zT2dist(zT2dist > 11,:);
    end

    % Now compute depths closest to T2B samples
    nmrDepths = z_dk;
    interpT2B = interp1(T2B_depth, T2B_peak, nmrDepths);

    % Match depths of permeability measurements from DPP data
    % with the depths in the new processed data -- need to find nearest
    % neighbors, exlude NMR points that do not have corresponding k
    % measurement. 

    for i = 1:length(Dk)
       [errorz(i), ind(i)] = min(abs(z_dk(i) - z)); 
    end


    %% Set rest of variables
%     phi = dparam(ind,2); % totalf
%     ClayH20 = dparam(ind,3); % clayf
%     CapH20 = dparam(ind,4); % capf
%     FreeH20 = dparam(ind,5); % freef
%     T2ML = dparam(ind,6); %mlT2
%     k_SDR_VC = dparam(ind,7); %Ksdr 
%     k_TC = dparam(ind,8); %Ktc 
%     k_SOE = dparam(ind,9); % Ksoe
%     Tsdr = dparam(ind,10); % Tsdr
%     Ttc = dparam(ind,11); % Ttc
%     Tsoe = dparam(ind,12); % Tsoe
%     soe = dparam(ind,13); % soe
%     noise = dparam(ind,14); %noise
%     %soe_3s = dparam(ind,9); 
    %soe_twm = dparam(ind,10); 
    %soe_twm_3s = dparam(ind,11); 

    z = z(ind); 

    % Estimate T-C parameters at DPP K intervals
    % Specify cutoff between FFI and BVI in ms
    %cutoff = 33*10^-3;

    [val, cutoffBin] = min(abs(T2linbins - cutoffs(kk)));
    zStep = 20;

    boundT2dist(1:cutoffBin) = T2dist(zStep,1:cutoffBin);
    boundT2dist(cutoffBin+1:100) = 0;

    freeT2dist(1:cutoffBin) = 0;
    freeT2dist(cutoffBin+1:100) = T2dist(zStep,cutoffBin+1:100);

    figure(1)
    hold on
    box on
    grid on
    plot(T2linbins, T2dist(zStep,:))
    plot([cutoffs(kk),cutoffs(kk)], [0,0.02],'LineWidth',2)
    set(gca, 'xscale','log')
    xlabel('T_2 (s)')
    ylabel('Relative Amplitude')
    set(gca,'FontSize',14)

    % plot(T2linbins,boundT2dist)
    % plot(T2linbins,freeT2dist)

    filtT2dist = T2dist(ind,:);

    % Signal amplitude is calibrated to sum to porosity
    % sum(boundT2dist) + sum(freeT2dist) = phi

    for k = 1:length(ind)

        boundT2dist(1:cutoffBin) = filtT2dist(k,1:cutoffBin);
        boundT2dist(cutoffBin+1:100) = 0;

        freeT2dist(1:cutoffBin) = 0;
        freeT2dist(cutoffBin+1:100) = filtT2dist(k,cutoffBin+1:100);

        BVI(k) = sum(boundT2dist);
        FFI(k) = sum(freeT2dist);     

        if BVI(k) == 0
            BVI(k) = 0.0001;
        end

    end
    
    indexQuotient = (FFI./BVI)';
    
    SDR_Kest = SDR_K(SDR_b(kk),SDR_m(kk),SDR_n(kk),phi,T2ML);
    %TC_Kest = TC_K(TC_c(kk),TC_m(kk),TC_n(kk),phi,indexQuotient);
    SOE_Kest = SOE_K(SOE_b(kk),SOE_n(kk),SumEch);
    Seevers_Kest = Seevers_K(Seevers_b(kk),Seevers_m(kk),Seevers_n(kk),...
        T2ML,interpT2B,phi);
    
    KGM_lKest = KGM_lK(KGM_rho(kk),KGM_tau(kk),KGM_m(kk),log10(phi),T2ML);
    KGM_Kest = 10.^KGM_lKest;
           
    SDR_error = computeError(Dk, SDR_Kest);
    %TC_error = computeError(Dk, TC_Kest);
    SOE_error = computeError(Dk, SOE_Kest);
    Seevers_error = computeError(Dk, Seevers_Kest);
    KGM_error = computeError(Dk, KGM_Kest);
    
    k_estimates = [SDR_Kest SOE_Kest Seevers_Kest KGM_Kest];
    k_names = {'DPP','SDR','TC 33ms','SOE','Seevers','KGM'};
    k_sym = {'+','o','*','x','^','v'};
    k_error = [SDR_error SOE_error Seevers_error KGM_error];
    
    errorEstimates(kk, :) = k_error;
    
    plotKestKdpp(K,k_estimates,k_names,k_sym)
    title(baseName)
    
    plotKwithDepth(K,z,[zT2dist T2dist],T2logbins,k_estimates,k_names,k_sym)
    title(baseName)

    clearvars z ind indexQuotient k_estimates phi T2ML FFI BVI boundT2dist freeT2dist

    
    
end