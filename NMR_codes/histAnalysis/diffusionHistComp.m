% Compare fast vs slow diffusion regime
% Plot K diff Histogram
sites = [{'Site1-WellG6'} {'Site1-WellG5'} {'Site2-WellPN1'} {'Site2-WellPN2'}];

sites_Maurer = {'dpnmr_larned_east','dpnmr_larned_lwph','dpnmr_larned_west',...
   'dpnmrA11','dpnmrA12','dpnmrC1S','dpnmrC1SE','dpnmrC1SW',...
   'dpnmr_leque_east','dpnmr_leque_west'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SDR_K = @(b,m,n,phi,T2ML) b.*(phi.^m).*(T2ML).^n;
Seevers_K = @(b,m,n,T2ML,T2B,phi) b.*(phi).^m.*((T2ML.^(-1) - T2B.^(-1)).^(-1)).^n;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T2B = 2.2293;

load SDR_bestFit_1101_m1_n2.mat
sites = siteList;
SDR_b_fast = squeeze(totalbMatrix(1,1,:))';
SDR_n_fast = n;
SDR_n_fast = repmat(n,length(siteList),1);
SDR_m_fast = squeeze(totalmMatrix(1,1,:))';

load Seevers_bestFit_1101_m1_n2_T2BAvg.mat
sites = siteList;
Seevers_b_fast = squeeze(totalbMatrix(1,1,:))';
Seevers_n_fast = n;
Seevers_n_fast = repmat(n,length(siteList),1);
Seevers_m_fast = squeeze(totalmMatrix(1,1,:))';

load SDR_bestFit_1101_m1_n1.mat
sites = siteList;
SDR_b_slow = squeeze(totalbMatrix(1,1,:))';
SDR_n_slow = n;
SDR_n_slow = repmat(n,length(siteList),1);
SDR_m_slow = squeeze(totalmMatrix(1,1,:))';

load Seevers_bestFit_1101_m1_n2_T2BAvg.mat
sites = siteList;
Seevers_b_slow = squeeze(totalbMatrix(1,1,:))';
Seevers_n_slow = n;
Seevers_n_slow = repmat(n,length(siteList),1);
Seevers_m_slow = squeeze(totalmMatrix(1,1,:))';


for kk = 1:length(sites)

    
    [T2dist, T2logbins, nmrName] = loadRawNMRdata(sites{kk});
   
    T2dist_matrix{kk} = T2dist;
    T2logbins_matrix{kk} = T2logbins;
    
    [d, K, T2ML, phi, z, SumEch, logK, logT2ML, logPhi, SumEch_3s, SumEch_twm, ...
    SumEch_twm_3s] = loadnmrdata2(nmrName);
    
    Kall{kk} = K;
    Phiall{kk} = phi;
    T2MLall{kk} = T2ML;
    SumEchAll{kk} = SumEch;
    
    SDR_Kest_fast{kk} = SDR_K(SDR_b_fast(kk),SDR_m_fast(kk),SDR_n_fast(kk),phi,T2ML);
    
    Seevers_Kest_fast{kk} = Seevers_K(Seevers_b_fast(kk),Seevers_m_fast(kk),Seevers_n_fast(kk),...
        T2ML,T2B,phi);
    
    SDR_Kest_slow{kk} = SDR_K(SDR_b_slow(kk),SDR_m_slow(kk),SDR_n_slow(kk),phi,T2ML);
    
    Seevers_Kest_slow{kk} = Seevers_K(Seevers_b_slow(kk),Seevers_m_slow(kk),Seevers_n_slow(kk),...
        T2ML,T2B,phi);
    
    [SDR_sign_fast{kk} SDR_diffFactor_fast{kk}] = estimateKdiffFactor_withSign(K,SDR_Kest_fast{kk},1);
    [Seevers_sign_fast{kk} Seevers_diffFactor_fast{kk}] = estimateKdiffFactor_withSign(K,Seevers_Kest_fast{kk},1);
    
    
    [SDR_sign_slow{kk} SDR_diffFactor_slow{kk}] = estimateKdiffFactor_withSign(K,SDR_Kest_slow{kk},1);
    [Seevers_sign_slow{kk} Seevers_diffFactor_slow{kk}] = estimateKdiffFactor_withSign(K,Seevers_Kest_slow{kk},1);
    
end

T2B = 2.31;

load SDR_maurer_bestFit_1101_m1_n2.mat

sites = siteList;
SDR_b_fast = squeeze(totalbMatrix(1,1,:))';
SDR_n_fast = n;
SDR_n_fast = repmat(n,length(siteList),1);
SDR_m_fast = squeeze(totalmMatrix(1,1,:))';

load Seevers_maurer_bestFit_1101_m1_n2_T2B2.31.mat
sites = siteList;
Seevers_b_fast = squeeze(totalbMatrix(1,1,:))';
Seevers_n_fast = n;
Seevers_n_fast = repmat(n,length(siteList),1);
Seevers_m_fast = squeeze(totalmMatrix(1,1,:))';

load SDR_maurer_bestFit_1101_m1_n1.mat

sites = siteList;
SDR_b_slow = squeeze(totalbMatrix(1,1,:))';
SDR_n_slow = n;
SDR_n_slow = repmat(n,length(siteList),1);
SDR_m_slow = squeeze(totalmMatrix(1,1,:))';

load Seevers_maurer_bestFit_1101_m1_n1_T2B2.31.mat
sites = siteList;
Seevers_b_slow = squeeze(totalbMatrix(1,1,:))';
Seevers_n_slow = n;
Seevers_n_slow = repmat(n,length(siteList),1);
Seevers_m_slow = squeeze(totalmMatrix(1,1,:))';

for kk = 1:length(sites_Maurer)
    [T2dist, T2logbins, nmrName] = loadRawNMRdata(sites_Maurer{kk});
   
    T2dist_matrix_Maurer{kk} = T2dist;
    T2logbins_matrix_Maurer{kk} = T2logbins;
    
    [d, K_Maurer, T2ML_Maurer, phi_Maurer, z, SumEch_Maurer, logK, logT2ML, logPhi, SumEch_3s, SumEch_twm, ...
    SumEch_twm_3s] = loadnmrdata2(nmrName);
    
    Kall_Maurer{kk} = K_Maurer;
    Phiall_Maurer{kk} = phi_Maurer;
    T2MLall_Maurer{kk} = T2ML_Maurer;
    
    SumEchAll_Maurer{kk} = SumEch_Maurer;
    
    SDR_Kest_Maurer_fast{kk} = SDR_K(SDR_b_fast(kk),SDR_m_fast(kk),SDR_n_fast(kk),phi_Maurer,T2ML_Maurer);
    Seevers_Kest_Maurer_fast{kk} = Seevers_K(Seevers_b_fast(kk),Seevers_m_fast(kk),Seevers_n_fast(kk),...
    T2ML_Maurer,T2B,phi_Maurer);

    SDR_Kest_Maurer_slow{kk} = SDR_K(SDR_b_slow(kk),SDR_m_slow(kk),SDR_n_slow(kk),phi_Maurer,T2ML_Maurer);
    Seevers_Kest_Maurer_slow{kk} = Seevers_K(Seevers_b_slow(kk),Seevers_m_slow(kk),Seevers_n_slow(kk),...
        T2ML_Maurer,T2B,phi_Maurer);

    [SDR_sign_Maurer_fast{kk} SDR_diffFactor_Maurer_fast{kk}] = estimateKdiffFactor_withSign(K_Maurer,SDR_Kest_Maurer_fast{kk},1);
    [Seevers_sign_Maurer_fast{kk} Seevers_diffFactor_Maurer_fast{kk}] = estimateKdiffFactor_withSign(K_Maurer,Seevers_Kest_Maurer_fast{kk},1);
    
    [SDR_sign_Maurer_slow{kk} SDR_diffFactor_Maurer_slow{kk}] = estimateKdiffFactor_withSign(K_Maurer,SDR_Kest_Maurer_slow{kk},1);
    [Seevers_sign_Maurer_slow{kk} Seevers_diffFactor_Maurer_slow{kk}] = estimateKdiffFactor_withSign(K_Maurer,Seevers_Kest_Maurer_slow{kk},1);
        
end

Kall_vec = vertcat(Kall{:});
Kall_Maurer_vec = vertcat(Kall_Maurer{:});
T2MLall_vec = vertcat(T2MLall{:});
T2MLall_Maurer_vec = vertcat(T2MLall_Maurer{:});

SDR_diffFactor_all_fast = vertcat(SDR_diffFactor_fast{:});
Seevers_diffFactor_all_fast = vertcat(Seevers_diffFactor_fast{:});
SDR_diffFactor_all_slow = vertcat(SDR_diffFactor_slow{:});
Seevers_diffFactor_all_slow = vertcat(Seevers_diffFactor_slow{:});

SDR_diffFactor_all_Maurer_fast = vertcat(SDR_diffFactor_Maurer_fast{:});
Seevers_diffFactor_all_Maurer_fast = vertcat(Seevers_diffFactor_Maurer_fast{:});
SDR_diffFactor_all_Maurer_slow = vertcat(SDR_diffFactor_Maurer_slow{:});
Seevers_diffFactor_all_Maurer_slow = vertcat(Seevers_diffFactor_Maurer_slow{:});

SDR_sign_all_fast = vertcat(SDR_sign_fast{:});
Seevers_sign_all_fast = vertcat(Seevers_sign_fast{:});
SDR_sign_all_slow = vertcat(SDR_sign_slow{:});
Seevers_sign_all_slow = vertcat(Seevers_sign_slow{:});

SDR_sign_all_Maurer_fast = vertcat(SDR_sign_Maurer_fast{:});
Seevers_sign_all_Maurer_fast = vertcat(Seevers_sign_Maurer_fast{:});
SDR_sign_all_Maurer_slow = vertcat(SDR_sign_Maurer_slow{:});
Seevers_sign_all_Maurer_slow = vertcat(Seevers_sign_Maurer_slow{:});

SDR_diffFactor_total_fast = [SDR_diffFactor_all_fast; SDR_diffFactor_all_Maurer_fast];
Seevers_diffFactor_total_fast = [Seevers_diffFactor_all_fast; Seevers_diffFactor_all_Maurer_fast];
SDR_diffFactor_total_slow = [SDR_diffFactor_all_slow; SDR_diffFactor_all_Maurer_slow];
Seevers_diffFactor_total_slow = [Seevers_diffFactor_all_slow; Seevers_diffFactor_all_Maurer_slow];

SDR_sign_total_fast = [SDR_sign_all_fast; SDR_sign_all_Maurer_fast];
Seevers_sign_total_fast = [Seevers_sign_all_fast; Seevers_sign_all_Maurer_fast];
SDR_sign_total_slow = [SDR_sign_all_slow; SDR_sign_all_Maurer_slow];
Seevers_sign_total_slow = [Seevers_sign_all_slow; Seevers_sign_all_Maurer_slow];

% SDR_diffFactor_total = [SDR_diffFactor_all];
% Seevers_diffFactor_total = [Seevers_diffFactor_all];
% SOE_diffFactor_total = [SOE_diffFactor_all];
% KGM_diffFactor_total = [KGM_diffFactor_all];
% 
% SDR_sign_total = [SDR_sign_all];
% Seevers_sign_total = [Seevers_sign_all];
% SOE_sign_total = [SOE_sign_all];
% KGM_sign_total = [KGM_sign_all];


% Manually assign diff factor for a factor above 20
SDR_diffFactor_total_fast(SDR_diffFactor_total_fast >= 100) = 18;
Seevers_diffFactor_total_fast(Seevers_diffFactor_total_fast >= 100) = 18;

SDR_diffFactor_total_slow(SDR_diffFactor_total_slow >= 100) = 18;
Seevers_diffFactor_total_slow(Seevers_diffFactor_total_slow >= 100) = 18;
%% 
% Plot histogram of Kdiff Factor
%edges = [1 1.5 2 2.5 3 4 5 10 20 40 80];
edgesLog = [-2 -1 -0.8 -0.6 -0.4 -0.2 0 0.2 0.4 0.6 0.8 1 2];
edges = 10.^edgesLog;

figure(4)

subplot(2,2,1)
hold on
grid on
box on
grid minor

title('SDR n = 2')

histDataSDR_fast = SDR_sign_total_fast.*log10(SDR_diffFactor_total_fast);

%histogram(SDR_diffFactor_all, edges)
histogram(histDataSDR_fast,edgesLog)
ylim([0 50])
%ylim([0 30])

xticklabels({'-100','-10','0','10','100'})
xlabel('K Difference Factor')
ylabel('Counts')
set(gca,'FontSize',14)

subplot(2,2,2)
hold on
grid on
box on
grid minor

title('SDR n = 1')
histDataSDR_slow = SDR_sign_total_slow.*log10(SDR_diffFactor_total_slow);

%histogram(Seevers_diffFactor_all, edges)
histogram(histDataSDR_slow,edgesLog)
ylim([0 50])
%ylim([0 30])

xticklabels({'-100','-10','0','10','100'})
xlabel('K Difference Factor')
ylabel('Counts')
set(gca,'FontSize',14)


subplot(2,2,3)

hold on
grid on
box on
grid minor

title('Seevers n = 2')
histDataSeevers_fast = Seevers_sign_total_fast.*log10(Seevers_diffFactor_total_fast);


%histogram(SOE_diffFactor_all, edges)
histogram(histDataSeevers_fast,edgesLog)
ylim([0 50])
%ylim([0 30])

xticklabels({'-100','-10','0','10','100'})
xlabel('K Difference Factor')
ylabel('Counts')
set(gca,'FontSize',14)


subplot(2,2,4)
hold on
grid on
box on
grid minor

title('Seevers n = 1')
histDataSeevers_slow = Seevers_sign_total_slow.*log10(Seevers_diffFactor_total_slow);


%histogram(KGM_diffFactor_all, edges)
histogram(histDataSeevers_slow,edgesLog)
ylim([0 50])
%ylim([0 25])

xticklabels({'-100','-10','0','10','100'})
xlabel('K Difference Factor')
ylabel('Counts')
set(gca,'FontSize',14)

