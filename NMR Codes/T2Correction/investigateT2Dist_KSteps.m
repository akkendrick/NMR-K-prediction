% Attempt to fit the data with two models of K

sites = [{'Site1-WellG6'} {'Site1-WellG5'} {'Site2-WellPN1'} {'Site2-WellPN2'}];

sites_Maurer = {'dpnmr_larned_east','dpnmr_larned_lwph','dpnmr_larned_west',...
   'dpnmrA11','dpnmrA12','dpnmrC1S','dpnmrC1SE','dpnmrC1SW',...
   'dpnmr_leque_east','dpnmr_leque_west'};

%sites_Maurer = {'dat

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
end

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

worstPoints_Maurer_Ind = find(KdiffFactor_Maurer > 10);
worstPoints_Ind = find(KdiffFactor > 10);

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
% Try including KGM
Temp = 20;  % temperature in degress C 
rho = @(Tt) 1000*(1 - ((Tt+288.94)./(508929*(Tt+68.12))).*(Tt-3.98).^2); % kg/m^3
eta = @(Tt) 0.0013 - 1.7e-5*Tt;         % Pa -s
Tb = @(Tt) 3.3 + 0.044*(Tt - 35);       % seconds
D = @(Tt) (1.0413 + 0.039828*Tt + 0.00040318*Tt.^2).*1e-9;  % m^2/s 
g = 9.8;    %m/s^2
tort = 1/(1.5^2); 
t1 = (rho(Temp)*g)/(8*eta(Temp)); % 

num2 = @(T2) 4*D(Temp)*Tb(Temp)*T2;
denom2 = @(T2) Tb(Temp) - T2; 
    
f12 = @(rho) (D(Temp)./rho);  
SQterm = @(rho,T2) sqrt(f12(rho).^2 + (num2(T2)./denom2(T2))); 

KGM_lK = @(rho,tau,m,lphi,T2) log10(1/tau^2) + log10(t1) + m*lphi + 2*log10(SQterm(rho,T2)-f12(rho)); 

% KGM (has been checked, looks good)
KGM_tau = [1 1.7378 2.7227 3.4674];
KGM_rho = [2.0606e-05 5.5335e-05 0.0077 1];
KGM_m = [0 0 0 0];

% KGM (has been checked, looks good)
KGM_tau_Maurer = [1.6788 1.9498 1 1 1.1482 1.5311 1.4289 1.4125 1 1.3335];
KGM_rho_Maurer = [0.0100 4.2170e-4 5.2723e-05 9.5499e-05 1.4191e-04 100 100 100 1.6141e-04 1.4191e-04];
KGM_m_Maurer = [0 0 0 0 0 0 0 0 0 0];

for kk = 1:length(sites)
    KGM_logK{kk} = KGM_lK(KGM_rho(kk), KGM_tau(kk), KGM_m(kk), Phiall{kk}, T2MLall{kk});
    KGM_K{kk} = 10.^KGM_logK{kk};
end

for kk = 1:length(sites_Maurer)
    KGM_logK_Maurer{kk} = KGM_lK(KGM_rho_Maurer(kk), KGM_tau_Maurer(kk), KGM_m_Maurer(kk), Phiall_Maurer{kk}, T2MLall_Maurer{kk});
    KGM_K_Maurer{kk} = 10.^KGM_logK_Maurer{kk};
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1)
box on
grid on
hold on

scatter(T2MLall_vec.^2, Kall_vec,40,'filled')
scatter(T2MLall_Maurer_vec.^2, Kall_Maurer_vec,40,'filled')

%scatter(T2SeeversAll_vec.^2, Kall_vec,40,'filled')
%scatter(T2SeeversAll_Maurer_vec.^2, Kall_Maurer_vec,40,'filled')

scatter(T2MLall_vec.^2, Kmodel_SDR_fast,20,'filled')
scatter(T2MLall_Maurer_vec.^2, Kmodel_Maurer_SDR_fast,20,'filled')

%scatter(T2MLall_vec.^2, Kmodel_Seevers,20,'filled')
%scatter(T2MLall_Maurer_vec.^2, Kmodel_Maurer_Seevers,20,'filled')

scatter((T2MLall_vec(worstPoints_Ind)).^2,Kall_vec(worstPoints_Ind),60,'filled','k')
scatter((T2MLall_Maurer_vec(worstPoints_Maurer_Seevers_Ind)).^2,Kall_Maurer_vec(worstPoints_Maurer_Seevers_Ind),60,'filled','k')

set(gca,'XScale','log')
set(gca,'YScale','log')
set(gca,'FontSize',14)
ylabel('Hydraulic Conductivity (m/s)')
xlabel('T2ML^2')

legend('Wisconsin Data','Kansas + Washington Data','Wisc SDR Model (n=2,m=0)','K+W SDR Model (n=2,m=0)','K Difference Factor > 10')

xlim([10^-4, 0.3])
ylim([10^-7, 10^-2])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(3)
box on
grid on
hold on

scatter(T2SeeversAll_vec.^2, Kall_vec,40,'filled')
scatter(T2SeeversAll_Maurer_vec.^2, Kall_Maurer_vec,40,'filled')

%scatter(T2SeeversAll_vec.^2, Kall_vec,40,'filled')
%scatter(T2SeeversAll_Maurer_vec.^2, Kall_Maurer_vec,40,'filled')

%scatter(T2MLall_vec.^2, Kmodel_SDR_fast,20,'filled')
%scatter(T2MLall_Maurer_vec.^2, Kmodel_Maurer_SDR_fast,20,'filled')

scatter(T2SeeversAll_vec.^2, Kmodel_Seevers,20,'filled')
scatter(T2SeeversAll_Maurer_vec.^2, Kmodel_Maurer_Seevers,20,'filled')

scatter((T2SeeversAll_vec(worstPoints_Seevers_Ind)).^2,Kall_vec(worstPoints_Seevers_Ind),60,'filled','k')
scatter((T2SeeversAll_Maurer_vec(worstPoints_Maurer_Ind)).^2,Kall_Maurer_vec(worstPoints_Maurer_Ind),60,'filled','k')

set(gca,'XScale','log')
set(gca,'YScale','log')
set(gca,'FontSize',14)
ylabel('Hydraulic Conductivity (m/s)')
xlabel('T2ML^2')

legend('Wisconsin Data','Kansas + Washington Data','Wisc SDR Model (n=2,m=0)','K+W SDR Model (n=2,m=0)','K Difference Factor > 10')

xlim([10^-4, 0.3])
ylim([10^-7, 10^-2])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Try to investigate T2 distribution of points with K diff > 10

%T2 dist with K diff > 10
%T2dist_filt
%T2logbins_filt

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% Investigate KGM
KGM_KAll = vertcat(KGM_K{:});
KGM_KAll_Maurer = vertcat(KGM_K_Maurer{:});

KdiffFactor_Maurer = estimateKdiffFactor(Kall_Maurer_vec,KGM_KAll_Maurer,1);
KdiffFactor = estimateKdiffFactor(Kall_vec,KGM_KAll,1);

worstPoints_Maurer_Ind = find(KdiffFactor_Maurer > 10);
worstPoints_Ind = find(KdiffFactor > 10);

figure(4)
box on
grid on
hold on

scatter(T2MLall_vec.^2, Kall_vec,40,'filled')
scatter(T2MLall_Maurer_vec.^2, Kall_Maurer_vec,40,'filled')

%scatter(T2SeeversAll_vec.^2, Kall_vec,40,'filled')
%scatter(T2SeeversAll_Maurer_vec.^2, Kall_Maurer_vec,40,'filled')

%scatter(T2MLall_vec.^2, Kmodel_SDR_fast,20,'filled')
%scatter(T2MLall_Maurer_vec.^2, Kmodel_Maurer_SDR_fast,20,'filled')

scatter(T2MLall_vec.^2, KGM_KAll,20,'filled')
scatter(T2MLall_Maurer_vec.^2, KGM_KAll_Maurer,20,'filled')

scatter((T2MLall_vec(worstPoints_Seevers_Ind)).^2,Kall_vec(worstPoints_Seevers_Ind),60,'filled','k')
scatter((T2MLall_Maurer_vec(worstPoints_Maurer_Ind)).^2,Kall_Maurer_vec(worstPoints_Maurer_Ind),60,'filled','k')

set(gca,'XScale','log')
set(gca,'YScale','log')
set(gca,'FontSize',14)
ylabel('Hydraulic Conductivity (m/s)')
xlabel('T2ML^2')

legend('Wisconsin Data','Kansas + Washington Data','Wisc SDR Model (n=2,m=0)','K+W SDR Model (n=2,m=0)','K Difference Factor > 10')

xlim([10^-4, 0.3])
ylim([10^-7, 10^-2])

