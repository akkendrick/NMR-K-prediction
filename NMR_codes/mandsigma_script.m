%% Script to estimate m and sigma at each site

% use each individual well
all_names = {'A1', 'C1', 'dpnmr_leque_east', 'dpnmr_leque_west', ...
  'dpnmr_larned_east', 'dpnmr_larned_west', 'dpnmr_larned_lwph'}; 

% all_names = {'G6_W2_tr5_20x_16p75_up_F_wRIN_wRFI_reg50_Va1', ...
%     'G5_W1_tr5_20x_16p5_up_F1n2_wRIN_wRFI_Reg50_Va1', ...
%     'Pl_W1_Tr5_20x_MPp75aLS_F1n2_wRIN_wRFI_Reg50_Va1', ...
%     'W2_Tr5_20x_MPp75aLS_Reg50_wRIN_wRFI_Va1'};

% use each site
% all_names = {'dpnmr_larned_all', 'dpnmr_leque_all', 'gems_all'}; 

% use all data together
% all_names = {'all_data'}; 
figureson = 0; 

% [b_mcmc, n_mcmc, m_mcmc, sig_mcmc] = deal(zeros(
for k = 1:length(all_names)
    name = all_names{k}
    [d, Dk, T2ML, phi, z, SumEch, kk, lt, lp, SumEch_3s, SumEch_twm, ...
        SumEch_twm_3s] = loadnmrdata2(name); 

    %% MCMC for solution to various parameters
    Niter= 1e6; 
    stepsize = 0.8; 

    % MCMC over all parameters (T2B, b, n, m, data error (sigma_mcmc)
    [paramhats, ~, ~, ~, accept_rat] = mcmc_nmr_full(Dk, T2ML, phi, ...
        z, Niter, stepsize, figureson);
    
    b_mcmc(:,k) = paramhats(2,:); 
    n_mcmc(:,k) = paramhats(3,:);
    m_mcmc(:,k) = paramhats(4,:);
    sig_mcmc(:,k) = paramhats(5,:);
    
    % %%%%%%%%%%%%%% MCMC for b and data error only. n and m are fixed 
    n = 2; 
    m = 4;

    [b,sig,~, accept_rat] = mcmc_nmr_bsig(Dk, T2ML, phi, z, m, n, ...
        Niter, stepsize, figureson);
    
    b_fixed(:,k) = b;
    sig_fixed(:,k) = sig;
    
    
    
end

%% Compute b statistics from fixed n and m data
logb_mean = mean(b_fixed);
b_mean = 10.^logb_mean


%% Plotting
dnames = {'A1', 'C1', 'Leque East', 'Leque West','Larned East', ...
    'Larned West', 'Larned C'}; 

%dnames = {'G6','G5','PN1','PN2'};

figure; 
subplot(322)
for k = 1:length(all_names)
    histogram(n_mcmc(:,k), 'Normalization', 'pdf')
    hold on
end
xlim([0,4])
xlabel('n', 'FontSize', 14)
ylabel('Density', 'FontSize', 14)
set(gca, 'FontSize', 12)

subplot(321)
for k = 1:length(all_names)
    histogram(m_mcmc(:,k), 'Normalization', 'pdf')
    hold on
end
xlabel('m', 'FontSize', 14)
ylabel('Density', 'FontSize', 14)
xlim([-2.1, 5])
set(gca, 'FontSize', 12)

subplot(323)
for k = 1:length(all_names)
    histogram(sig_mcmc(:,k), 'Normalization', 'pdf')
    hold on
end
%xlim([0.1, 1.8])
xlim([0.1, 2])

xlabel('\sigma', 'FontSize', 14)
ylabel('Density', 'FontSize', 14)
set(gca, 'FontSize', 12)
%legend(dnames, 'FontSize', 14)

subplot(324)
for k = 1:length(all_names)
    histogram(b_mcmc(:,k), 'Normalization', 'pdf')
    hold on
end
xlim([-6, 4])
xlabel('b', 'FontSize', 14)
ylabel('Density', 'FontSize', 14)
set(gca, 'FontSize', 12)

subplot(325)
for k = 1:length(all_names)
    histogram(sig_fixed(:,k), 'Normalization', 'pdf')
    hold on
end
%xlim([0.1, 1.8])
xlim([0.1, 2])

xlabel('\sigma Fixed', 'FontSize', 14)
ylabel('Density', 'FontSize', 14)
set(gca, 'FontSize', 12)

subplot(326)
for k = 1:length(all_names)
    histogram(b_fixed(:,k), 'Normalization', 'pdf')
    hold on
end
%xlim([-3.2, 0])
%xlim([-3.2, -1.5])
xlim([-6, 4])
xlabel('b fixed', 'FontSize', 14)
ylabel('Density', 'FontSize', 14)
set(gca, 'FontSize', 12)
legend(dnames, 'FontSize', 14)
