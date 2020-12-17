function [K,z,T2dist,T2logbins,kTC_best,bestFitMatrix,totalError,meanErrorFactor,medianErrorFactor,indexQuotient] = computeTCperm(baseName,n,m,cutoff,figureson)

%baseDir = '/Volumes/GoogleDrive/My Drive/USGS Project/USGS Data/';
baseDir = 'I:\My Drive\Stanford\USGS Project\Field Data\USGS Data\';
%baseDir = 'I:\My Drive\Stanford\USGS Project\Kansas_Wash_Data\';
%baseDir = '/Volumes/GoogleDrive/My Drive/USGS Project/Kansas_Wash_Data/';

if strcmp(baseName,'Site1-WellG5')
    site = 'Site1-WellG5';
    name = 'G5_W1_tr5_20x_16p5_up_F1n2_wRIN_wRFI_Reg50_Va1';
    nmrName = name;
elseif strcmp(baseName,'Site1-WellG6')
    site = 'Site1-WellG6';
    name = 'G6_W2_tr5_20x_16p75_up_F_wRIN_wRFI_reg50_Va1';
    nmrName = name;
elseif strcmp(baseName,'Site2-WellPN1')
    site = 'Site2-WellPN1';
    name = 'Pl_W1_Tr5_20x_MPp75aLS_F1n2_wRIN_wRFI_Reg50_Va1';
    nmrName = name;
elseif strcmp(baseName,'Site2-WellPN2')
    site = 'Site2-WellPN2';
    name = 'W2_Tr5_20x_MPp75aLS_Reg50_wRIN_wRFI_Va1';
    nmrName = name;
else 
    site = baseName;
end

% %in3 = [baseDir site '/' name '/' strcat(site,'_DPP_filt.txt')];
% in5 = [baseDir site '/' site '_T2_dist.txt'];
% in6 = [baseDir site '/' site '_T2_bins_log10s.txt'];
% 
% 
% %DPPdat = load(in3); 
% T2dist = load(in5);
% T2logbins = load(in6);
% % set needed variables

%Dk = olddat(:,4)*1.16e-5; 
%z_dk = olddat(:,1); 

[T2dist, T2logbins, SEdecayTime,SEdecayUniform,SEdecay,oneDVectors,...
    oneDVectorsUniform, nmrName] = loadAllRawNMRdata(site);

[d, K, T2ML, phi, z, SumEch, logK, logT2ML, logPhi, SumEch_3s, SumEch_twm, ...
SumEch_twm_3s] = loadnmrdata2(site); 

Dk = K;
z_dk = z;

origT2dist = T2dist;

zT2dist = T2dist(:,1);
T2dist = T2dist(:,2:end);


T2linbins = 10.^T2logbins;

% Match depths of permeability measurements from DPP data
% with the depths in the new processed data -- need to find nearest
% neighbors, exlude NMR points that do not have corresponding k
% measurement. 

errorz = zeros(1,length(Dk));
ind = zeros(1,length(Dk));

for i = 1:length(Dk)
   [errorz(i), ind(i)] = min(abs(z_dk(i) - zT2dist)); 
end


%% Set rest of variables
% phi = dparam(ind,2); % totalf
% ClayH20 = dparam(ind,3); % clayf
% CapH20 = dparam(ind,4); % capf
% FreeH20 = dparam(ind,5); % freef
% T2ML = dparam(ind,6); %mlT2
% k_SDR_VC = dparam(ind,7); %Ksdr 
% k_TC = dparam(ind,8); %Ktc 
% k_SOE = dparam(ind,9); % Ksoe
% Tsdr = dparam(ind,10); % Tsdr
% Ttc = dparam(ind,11); % Ttc
% Tsoe = dparam(ind,12); % Tsoe
% soe = dparam(ind,13); % soe
% noise = dparam(ind,14); %noise
%soe_3s = dparam(ind,9); 
%soe_twm = dparam(ind,10); 
%soe_twm_3s = dparam(ind,11); 

z = zT2dist(ind); 

% Estimate T-C parameters at DPP K intervals
% Specify cutoff between FFI and BVI in ms
%cutoff = 33*10^-3;

[~, cutoffBin] = min(abs(T2linbins - cutoff));
% zStep = 20;
% 
% boundT2dist(1:cutoffBin) = T2dist(zStep,1:cutoffBin);
% boundT2dist(cutoffBin+1:100) = 0;
% 
% freeT2dist(1:cutoffBin) = 0;
% freeT2dist(cutoffBin+1:100) = T2dist(zStep,cutoffBin+1:100);
% 


% plot(T2linbins,boundT2dist)
% plot(T2linbins,freeT2dist)

filtT2dist = T2dist(ind,:);

% Signal amplitude is calibrated to sum to porosity
% sum(boundT2dist) + sum(freeT2dist) = phi

BVI = zeros(1,length(ind));
FFI = zeros(1,length(ind));
nBins = length(T2logbins);
for k = 1:length(ind)
    
    boundT2dist(1:cutoffBin) = filtT2dist(k,1:cutoffBin);
    boundT2dist(cutoffBin+1:nBins) = 0;

    freeT2dist(1:cutoffBin) = 0;
    freeT2dist(cutoffBin+1:nBins) = filtT2dist(k,cutoffBin+1:nBins);
    
    BVI(k) = sum(boundT2dist);
    
    FFI(k) = sum(freeT2dist);     
    
    if BVI(k) == 0
        BVI(k) = (randi(100)/1000)^randi(100);
    end
    
    if FFI(k) == 0
        FFI(k) = randi(100)*randi(100);
    end
    
    if figureson == 1
        figure(1)
        hold on
        box on
        grid on

        plot(T2linbins, filtT2dist(k,:))
        plot(T2linbins, boundT2dist)
        plot(T2linbins, freeT2dist)

        plot([cutoff,cutoff], [0,0.02],'LineWidth',2)
        set(gca, 'xscale','log')
        xlabel('T_2 (s)')
        ylabel('Relative Amplitude')
        set(gca,'FontSize',14)
        
        close all
   end

end

% cutoff
% BVI
% median(BVI)
% mean(BVI)
% std(BVI)
% 
% FFI

lkTC = @(c,m,n,lPhi,logFrac) log10(c) + m*lPhi + n*(logFrac);

% logkTC = log10(c) + m*lPhi + n*log10(FFI/BVI)

% Now that we have estimated BVI and FFI via cutoff, use bootstrap to
% estimate empirical parameters
indexQuotient = log10(FFI./BVI)';

Nboot = 2000;
%
% n can vary
%[c_boot, n_boot, m_boot] = bootstrap_fun([indexQuotient, logPhi, logK], Nboot,n);

if isempty(n) && isempty(m)
    % m, n can vary
    [c_boot, n_boot, m_boot] = bootstrap_fun([indexQuotient, logPhi, logK], Nboot);
else
    % n and m fixed
    [c_boot, n_boot, m_boot] = bootstrap_fun([indexQuotient, logPhi, logK], Nboot, n, m);
end

%cs = log10(c_boot); 
%graph_correlations([cs,n_boot,m_boot], 3, {'log_{10}(c)', 'n','m'}, 1, 0)
%graph_correlations([cs,m_boot], 3, {'log_{10}(c)','m'}, 1, 0)
%graph_correlations([cs, n_boot], 3, {'log_{10}(c)', 'n'}, 1, 0)

median_c = median(c_boot);
median_n = median(n_boot);
median_m = median(m_boot);

bestFitMatrix(1,1) = median_c;
bestFitMatrix(2,1) = median_n;
bestFitMatrix(3,1) = median_m;

lkTC_best = lkTC(median_c,median_m,median_n,logPhi,indexQuotient);
kTC_best = 10.^lkTC_best;
kTC_best;

%lkTC_Dlubac = kTC(1.6*10^-5,2,2,logPhi,indexQuotient);
%kTC_Dlubac = 10.^lkTC_Dlubac;

%plotKwithDepth(K,z,origT2dist,T2logbins,[kTC_best, kTC_Dlubac],...
%    {'DPP','T-C 33ms','T-C Dlubac'},{'+','o'})

% if figureson == 1
%     plotKwithDepth(K,z,origT2dist,T2logbins,[kTC_best],...
%         {'DPP','T-C 33ms'},{'+'})
% end

totalError = computeError(K,kTC_best);
meanErrorFactor = mean(estimateKdiffFactor(K,kTC_best,1));
medianErrorFactor = median(estimateKdiffFactor(K,kTC_best,1));

%dlubacError = computeError(K,kTC_Dlubac);




end

