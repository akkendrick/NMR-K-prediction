% Range over pairs of m and n values
%close all
clear

load enso 

figureson = 0;


% Kansas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  siteList = [{'dpnmr_larned_east'} {'dpnmr_larned_lwph'} {'dpnmr_larned_west'}... 
%      {'dpnmrA11'} {'dpnmrA12'} {'dpnmrC1S'} {'dpnmrC1SE'} {'dpnmrC1SW'}];
% 
% refCutoff = [33 33 33 33 33 33 33 33];
% idealCutoff = [84 20 46 40 52 164 174 92];
% 
% SDR_b = [0.027 0.016 0.024 0.032 0.032 0.028 0.028 0.028];
% SDR_m = [0 0 0 0 0 0 0 0];
% SDR_n = [2 2 2 2 2 2 2 2];
  
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Washington
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% siteList = [{'dpnmr_leque_east'}];% {'dpnmr_leque_west'}];
% 
% refCutoff = [33 33]./1000;
% idealCutoff = [98 106]./1000;
% 
% SDR_b = [0.036 0.022];
% SDR_m = [0 0];
% SDR_n = [2 2];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Wisc
siteList = [{'Site1-WellG5'} {'Site1-WellG6'},...
    {'Site2-WellPN1'} {'Site2-WellPN2'} ];

refCutoff = [33, 33, 33, 33]./1000;
%refCutoff = [3, 3, 3, 3]./1000;

%idealCutoff = [400, 380, 70, 100]./1000;
%idealCutoff = [931.4406, 598.5475, 247.2, 247.2]/1000; %mean error factor
%idealCutoff = [322.2697, 173.5163, 247.1642, 207.0918]/1000; %sum(RMSE) 5/22/20
% idealCutoff = [102.0640, 352.0714, 189.5622, 173.5163]/1000; %sum(RMSE)
idealCutoff = [42.1,42.1,42.1,42.1]./1000; % from median min of T2 gap 11/18/20, m = 0, n = 2
TC_c = [1.02*10^-08, 2.09*10^-09, 9.39*10^-09,7.81*10^-09]; 

%idealCutoff = [1017.6, 1449.5, 2.7165, 322.2697]./1000;

%idealCutoff = [1060, 548, 226, 318]./1000;
%idealCutoff = [1000 1000 1000 1000]./1000;

SDR_b = [0.0024 0.0040 0.0075 0.0047];
SDR_m = [0 0 0 0];
SDR_n = [2 2 2 2];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SDR_K = @(b,m,n,phi,T2ML) b.*(phi.^m).*(T2ML).^n;

%load('optimalCutoffTable_n2_m1_RMSE_2000.mat')
%idealCutoff = [1730, 598.5475,1.9071,270.0206]./1000;
%idealCutoff = [1017.6, 1449.5, 2.7165, 322.2697]./1000;


for kk = 1:length(siteList)
   
    [T2dist,T2logbins,SEdecayTime,SEdecayUniform,SEdecay,oneDVectors,...
    oneDVectorsUniform, nmrName] = loadAllRawNMRdata(siteList{kk});

    [d, K, T2ML, phi, z, SumEch, log10K, log10T2, log10Porosity,...
        SumEch_3s, SumEch_twm, SumEch_twm_3s] = loadnmrdata2(siteList{kk});
    
    phiMatrix{kk} = phi;
    T2MLMatrix{kk} = T2ML;
    KMatrix{kk} = K;
    zMatrix{kk} = z;
    
    [K,z,T2dist,T2logbins,kTC_ref_best{kk},bestFitMatrix,totalError,...
        indexQuotient] = computeTCperm_nobtstp(siteList{kk},SDR_n(kk),SDR_m(kk),TC_c(kk),refCutoff(kk),figureson);
    
    [K_best,z_best,T2dist_best,T2logbins_best,kTC_ideal_best{kk},bestFitMatrix_best,totalError_best,...
        indexQuotient_best] = computeTCperm_nobtstp(siteList{kk},SDR_n(kk),SDR_m(kk),TC_c(kk),idealCutoff(kk),figureson);

    SDR_Kest{kk} = SDR_K(SDR_b(kk),SDR_m(kk),SDR_n(kk),phiMatrix{kk},T2MLMatrix{kk});  
    
    [SDR_sign{kk} SDRError{kk}] = estimateKdiffFactor_withSign(KMatrix{kk},SDR_Kest{kk},1);  
    
    [TC_refSign{kk} TC_refError{kk}] = estimateKdiffFactor_withSign(KMatrix{kk},kTC_ref_best{kk},1);  
    [TC_idealSign{kk} TC_idealError{kk}] = estimateKdiffFactor_withSign(KMatrix{kk},kTC_ideal_best{kk},1);  
    
%     meanSDRError(kk) = mean(SDRError{kk});
%     meanTC_refError(kk) = mean(TC_refError{kk});
%     meanTC_idealError(kk) = mean(TC_idealError{kk});    
    
    SDRError{kk}(SDRError{kk} > 100) = 100;
    TC_refError{kk}(TC_refError{kk} > 100) = 100;
    TC_idealError{kk}(TC_idealError{kk} > 100) = 100;

    SDRError{kk} = log10(SDRError{kk}) .* SDR_sign{kk};
    TC_refError{kk} = log10(TC_refError{kk}) .* TC_refSign{kk};
    TC_idealError{kk} =  log10(TC_idealError{kk}) .*  TC_idealSign{kk};

    TCIdealK{kk} = kTC_ideal_best{kk};
    Kall{kk} = K;
end


% figure(1)
% KmodelDiffHist_TC(SDRError{1}, TC_refError{1},TC_idealError{1})
% 
% figure(2)
% KmodelDiffHist_TC(SDRError{2}, TC_refError{2},TC_idealError{2})
% 
% figure(3)
% KmodelDiffHist_TC(SDRError{3}, TC_refError{3},TC_idealError{3})
% 
% figure(4)
% KmodelDiffHist_TC(SDRError{4}, TC_refError{4},TC_idealError{4})

SDR_error = vertcat(SDRError{:});

TC_ref_error = vertcat(TC_refError{:});
TC_ideal_error = vertcat(TC_idealError{:});

x0=10;
y0=10;
width=800;
height=1000;
set(gcf,'position',[x0,y0,width,height])

figure(1)
KmodelDiffHist_TC(SDR_error, TC_ref_error, TC_ideal_error)

sum(abs(SDR_error))
sum(abs(TC_ref_error))
sum(abs(TC_ideal_error))
