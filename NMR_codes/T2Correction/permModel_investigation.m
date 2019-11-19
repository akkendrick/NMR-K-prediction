% Attempt to fit the data with two models of K

sites = [{'Site1-WellG6'} {'Site1-WellG5'} {'Site2-WellPN1'} {'Site2-WellPN2'}];

sites_Maurer = {'dpnmr_larned_east','dpnmr_larned_lwph','dpnmr_larned_west',...
   'dpnmrA11','dpnmrA12','dpnmrC1S','dpnmrC1SE','dpnmrC1SW',...
   'dpnmr_leque_east','dpnmr_leque_west'};

for kk = 1:length(sites)
    [T2dist, T2logbins, nmrName] = loadRawNMRdata(sites{kk});
   
    [d, K, T2ML, phi, z, SumEch, logK, logT2ML, logPhi, SumEch_3s, SumEch_twm, ...
    SumEch_twm_3s] = loadnmrdata2(nmrName);
    
    Kall{kk} = K;
    Phiall{kk} = phi;
    T2MLall{kk} = T2ML;
    SumEchAll{kk} = SumEch;
end

for kk = 1:length(sites_Maurer)
    [T2dist, T2logbins, nmrName] = loadRawNMRdata(sites_Maurer{kk});
   
    [d, K_Maurer, T2ML_Maurer, phi_Maurer, z, SumEch_Maurer, logK, logT2ML, logPhi, SumEch_3s, SumEch_twm, ...
    SumEch_twm_3s] = loadnmrdata2(nmrName);
    
    Kall_Maurer{kk} = K_Maurer;
    Phiall_Maurer{kk} = phi_Maurer;
    T2MLall_Maurer{kk} = T2ML_Maurer;
    
    SumEchAll_Maurer{kk} = SumEch_Maurer;
        
end

Kall_vec = vertcat(Kall{:});
Kall_Maurer_vec = vertcat(Kall_Maurer{:});
T2MLall_vec = vertcat(T2MLall{:});
T2MLall_Maurer_vec = vertcat(T2MLall_Maurer{:});

bWisc = 0.0043;
bMaurer = 0.0278;
nModel = 2;

Kmodel_Maurer_SDR_fast = bMaurer.*(T2MLall_Maurer_vec).^(nModel);
Kmodel_SDR_fast = bWisc.*(T2MLall_vec).^(nModel);

KdiffFactor_Maurer = estimateKdiffFactor(Kall_Maurer_vec,Kmodel_Maurer_SDR_fast,1);
KdiffFactor = estimateKdiffFactor(Kall_vec,Kmodel_SDR_fast,1);

worstPoints_Maurer_Ind = find(KdiffFactor_Maurer > 8);
worstPoints_Ind = find(KdiffFactor > 8);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Try including Seevers
bWiscSeevers = 0.0036;
bMaurerSeevers = 0.0118;

T2B = 2.2;
T2SeeversAll_vec = ((1./T2MLall_vec - 1./T2B).^(-1));
T2SeeversAll_Maurer_vec = ((1./T2MLall_Maurer_vec - 1./T2B).^(-1));

Kmodel_Maurer_Seevers = bMaurerSeevers.*(T2SeeversAll_Maurer_vec).^(nModel);
Kmodel_Seevers = bWiscSeevers.*(T2SeeversAll_vec).^(nModel);

KdiffFactor_Maurer_Seevers = estimateKdiffFactor(Kall_Maurer_vec,Kmodel_Maurer_Seevers,1);
KdiffFactor_Seevers = estimateKdiffFactor(Kall_vec,Kmodel_Seevers,1);

worstPoints_Maurer_Seevers_Ind = find(KdiffFactor_Maurer_Seevers > 10);
worstPoints_Seevers_Ind = find(KdiffFactor_Seevers > 10);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1)
box on
grid on
hold on

scatter(T2MLall_vec.^2, Kall_vec,40,'filled')
scatter(T2MLall_Maurer_vec.^2, Kall_Maurer_vec,40,'filled')

%scatter(T2SeeversAll_vec.^2, Kall_vec,40,'filled')
%scatter(T2SeeversAll_Maurer_vec.^2, Kall_Maurer_vec,40,'filled')

%scatter(T2MLall_vec.^2, Kmodel_SDR_fast,20,'filled')
%scatter(T2MLall_Maurer_vec.^2, Kmodel_Maurer_SDR_fast,20,'filled')

scatter(T2MLall_vec.^2, Kmodel_Seevers,20,'filled')
scatter(T2MLall_Maurer_vec.^2, Kmodel_Maurer_Seevers,20,'filled')

scatter((T2MLall_vec(worstPoints_Ind)).^2,Kall_vec(worstPoints_Ind),60,'filled','k')
scatter((T2MLall_Maurer_vec(worstPoints_Maurer_Ind)).^2,Kall_Maurer_vec(worstPoints_Maurer_Ind),60,'filled','k')

set(gca,'XScale','log')
set(gca,'YScale','log')
set(gca,'FontSize',14)
ylabel('Hydraulic Conductivity (m/s)')
xlabel('T2ML^2')

legend('Wisconsin Data','Kansas + Washington Data','Wisc SDR Model (n=2,m=0)','K+W SDR Model (n=2,m=0)','K Difference Factor > 10')

xlim([10^-4, 0.3])
ylim([10^-7, 10^-2])

%%
% Investigate SOE

SumEchAll_vec = vertcat(SumEchAll{:});
SumEchAll_Maurer_vec = vertcat(SumEchAll_Maurer{:});

KSOEmodel = 0.0087.*vertcat(SumEchAll{:});
KSOEmodel_Maurer = 0.0043.*vertcat(SumEchAll_Maurer{:});

KdiffFactor_Maurer = estimateKdiffFactor(Kall_Maurer_vec,KSOEmodel_Maurer,1);
KdiffFactor = estimateKdiffFactor(Kall_vec,KSOEmodel,1);

worstPoints_Maurer_Ind = find(KdiffFactor_Maurer > 10);
worstPoints_Ind = find(KdiffFactor > 10);


figure(2)
box on
grid on
hold on

scatter(vertcat(SumEchAll{:}), Kall_vec,40,'filled')
scatter(vertcat(SumEchAll_Maurer{:}),Kall_Maurer_vec,40,'filled')

scatter(vertcat(SumEchAll{:}),KSOEmodel,20,'filled')
scatter(vertcat(SumEchAll_Maurer{:}), KSOEmodel_Maurer,20,'filled')

scatter((SumEchAll_vec(worstPoints_Ind)),Kall_vec(worstPoints_Ind),60,'filled','k')
scatter((SumEchAll_Maurer_vec(worstPoints_Maurer_Ind)),Kall_Maurer_vec(worstPoints_Maurer_Ind),60,'filled','k')


legend('Wisconsin Data','Kansas + Washington Data','Wisc SOE Model (n=1)','K+W SOE Model (n=1)','K Difference Factor > 10')

%scatter(T2MLall_vec.^2, Kmodel_SDR_fast,20,'filled')
%scatter(T2MLall_Maurer_vec.^2, Kmodel_Maurer_SDR_fast,20,'filled')

set(gca,'XScale','log')
set(gca,'YScale','log')
set(gca,'FontSize',14)
ylabel('Hydraulic Conductivity (m/s)')
xlabel('SOE')

xlim([0.001, 0.3])
ylim([10^-7, 10^-2])
