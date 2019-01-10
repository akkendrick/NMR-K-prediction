%% Script to estimate b directly

clear

%all_names = {'A1', 'C1', 'dpnmr_larned_east', 'dpnmr_larned_west', ...
%  'dpnmr_larned_lwph', 'dpnmr_leque_east', 'dpnmr_leque_west'}; 

all_names = { 'G5_W1_tr5_20x_16p5_up_F1n2_wRIN_wRFI_Reg50_Va1', ...
    'G6_W2_tr5_20x_16p75_up_F_wRIN_wRFI_reg50_Va1', ...
    'Pl_W1_Tr5_20x_MPp75aLS_F1n2_wRIN_wRFI_Reg50_Va1', ...
    'W2_Tr5_20x_MPp75aLS_Reg50_wRIN_wRFI_Va1'};


for k = 1:length(all_names)
   
    name = all_names{k}
    [d, K, T2ML, phi, z, SumEch, logK, logT2ML, logPhi, SumEch_3s, SumEch_twm, ...
        SumEch_twm_3s] = loadnmrdata2(name); 
    logSumEch = log10(SumEch); 


    figureson = 1; 

    %% Basic solving for b for fixed n, m
    % Given m and n, we can solve directly for b -> b = log(k) - m*log(phi) -
    % n*log(T2ML). This is the 'direct' method.

    C = @(m, n, lt, lp) m*lp + n*lt; 
    n = 2;
    m = 4; 
    bdir_n1 = logK - C(m, n, logT2ML, logPhi); 

    n = 2; 
    bdir_n2 = logK - C(m, n, logT2ML, logPhi); 

    if figureson == 1
        figure; 
        subplot(211)
        hist(bdir_n1, floor(sqrt(length(bdir_n1))))
        xlabel('b')
        ylabel('Frequency')
        title('b values estimated directly for fixed n, n = 1')
        subplot(212)
        hist(bdir_n2, floor(sqrt(length(bdir_n2))))
        xlabel('b')
        ylabel('Frequency')
        title('b values estimated directly for fixed n, n = 2')
    end
    
    logb_mean = mean(bdir_n2);
    b_mean(k) = 10.^logb_mean;

    logb_median = median(bdir_n2);
    b_median(k) = 10.^logb_median;
    
end

b_mean
b_median


