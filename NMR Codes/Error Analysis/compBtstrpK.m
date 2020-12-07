% Compare K estimates from each range of bootstrapped values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% K estimates
SDR_K = @(b,m,n,phi,T2ML) b.*(phi.^m).*(T2ML).^n;
SOE_K = @(b,n,SOE) b.*(SOE).^n;
Seevers_K = @(b,m,n,T2ML,T2B,phi) b.*(phi).^m.*((T2ML.^(-1) - T2B.^(-1)).^(-1)).^n;

Temp = 10.6;  % temperature in degress C 
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
siteList = [{'Site1-WellG5'}...
   {'Site1-WellG6'} {'Site2-WellPN1'} {'Site2-WellPN2'}];

for kk = 1:length(siteList)
    [T2dist{kk}, T2logbins{kk}, nmrName] = loadRawNMRdata(siteList{kk});
    depths{kk} = T2dist{kk}(:,1);
    T2dist{kk} = T2dist{kk}(:,2:end);
    
    [d{kk}, K{kk}, T2ML{kk}, phi{kk}, z{kk}, SumEch{kk},logK,...
        logT2,logPhi,sumech3,sumech_twm,sumechtwm3s] = loadnmrdata2(nmrName);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load KGM_btstrp_Wisc_1000_10.6DegC_m1.mat
m = 1;

nBoot = 1000;

for kk = 1:length(siteList)
  
   currentTau = bootKGM{1,kk}(:,1);
   currentRho = bootKGM{1,kk}(:,2);
   log10Rho = log10(currentRho);
   
   currentData = data{1,kk};
%    T2ML = 10.^(currentData(:,1));
   lphi = log10(currentData(:,2));
   logK = currentData(:,3);
   
   rhoBins =  linspace(-4, -1, 20); 
   tauBins = linspace(1,2.5, 40);

   medianlog10Rho = median(log10Rho(log10Rho < 2));
   medianTau = median(currentTau);
   
   KGM_Kmedian{kk} = 10.^(KGM_lK(10.^medianlog10Rho, medianTau,m,lphi,T2ML{kk}));
   
   for j = 1:nBoot
        KGM_Kestimates{j,kk} = 10.^(KGM_lK(currentRho(j), currentTau(j),m,lphi,T2ML{kk}));
        KGM_Kfactor{j,kk} = estimateKdiffFactor(K{kk},KGM_Kestimates{j,kk},1);
        KGM_meanKFactor{j,kk} = mean(KGM_Kfactor{j,kk});
   end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
load SDR_wisc_bboot_n2_m1.mat
nBoot = 2000;

n = 2;
m = 1;

SDR_K = @(b,m,n,phi,T2ML) b.*(phi.^m).*(T2ML).^n;

for kk = 1:length(siteList)

    currentb = b_boot_all{1,kk};
    medianb = median(b_boot_all{1,kk});
    
    SDR_Kmedian{kk} = SDR_K(medianb,m,n,phi{kk},T2ML{kk});
    
    for j = 1:nBoot
        SDR_Kestimates{j,kk} = SDR_K(currentb(j),m,n,phi{kk},T2ML{kk});
        SDR_Kfactor{j,kk} = estimateKdiffFactor(K{kk},SDR_Kestimates{j,kk},1);
        SDR_meanKFactor{j,kk} = mean(SDR_Kfactor{j,kk});
    end
    
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load Seevers_bestFit_1201_m1_n2_T2BAvg.mat
nBoot = 2000;

n = 2;
m = 1;

Seevers_K = @(b,m,n,T2ML,T2B,phi) b.*(phi).^m.*((T2ML.^(-1) - T2B.^(-1)).^(-1)).^n;

T2Bavg = 2.2293;

for kk = 1:length(siteList)

    currentb = b_boot_all{1,kk};
    medianb = median(b_boot_all{1,kk});
    
    Seevers_Kmedian{kk} = Seevers_K(medianb,m,n,T2ML{kk},T2Bavg,phi{kk});
    
    for j = 1:nBoot
        Seevers_Kestimates{j,kk} = Seevers_K(currentb(j),m,n,T2ML{kk},T2Bavg,phi{kk});
        Seevers_Kfactor{j,kk} = estimateKdiffFactor(K{kk},Seevers_Kestimates{j,kk},1);
        Seevers_meanKFactor{j,kk} = mean(Seevers_Kfactor{j,kk});
    end
    
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load SOE_n1_Wisc.mat
nBoot = 1000;

n = 1;
SOE_K = @(b,n,SOE) b.*(SOE).^n;

for kk = 1:length(siteList)
   
    currentb = b_boot_all{1,kk};
    medianb = median(b_boot_all{1,kk});
    
    SOE_Kmedian{kk} = SOE_K(medianb,n,SumEch{kk});
    
    for j = 1:nBoot
        SOE_Kestimates{j,kk} = SOE_K(currentb(j),n,SumEch{kk});
        SOE_Kfactor{j,kk} = estimateKdiffFactor(K{kk},SOE_Kestimates{j,kk},1);
        SOE_meanKFactor{j,kk} = mean(SOE_Kfactor{j,kk});
    end
    
end

%% Organize data by depth

for kk = 1:length(siteList)

    
    for j = 1:2000
        SDR_row(:,j) = SDR_Kestimates{j,kk};
    end

    logSDR_row = log10(SDR_row);

    maxlogK_row = max(logSDR_row,[],2);
    minlogK_row = min(logSDR_row,[],2);

    SDRKdiff = maxlogK_row - minlogK_row;
    meanSDRKdiff(kk) = mean(SDRKdiff)

    %%%%%

    for j = 1:2000
        Seevers_row(:,j) = Seevers_Kestimates{j,kk};
    end

    logSeevers_row = log10(Seevers_row);

    maxlogK_row = max(logSeevers_row,[],2);
    minlogK_row = min(logSeevers_row,[],2);

    SeeversKdiff = maxlogK_row - minlogK_row;
    meanSeeversKdiff(kk) = mean(SeeversKdiff)
    %%%%%

    for j = 1:1000
        SOE_row(:,j) = SOE_Kestimates{j,kk};
    end

    logSOE_row = log10(SOE_row);

    maxlogK_row = max(logSOE_row,[],2);
    minlogK_row = min(logSOE_row,[],2);

    SOEKdiff = maxlogK_row - minlogK_row;
    meanSOEKdiff(kk) = mean(SOEKdiff)

    %%%%%

    for j = 1:1000
        KGM_row(:,j) = KGM_Kestimates{j,kk};
    end

    logKGM_row = log10(KGM_row);

    maxlogK_row = max(logKGM_row,[],2);
    minlogK_row = min(logKGM_row,[],2);

    KGMKdiff = maxlogK_row - minlogK_row;

    meanKGMKdiff(kk) = mean(KGMKdiff)
    
    clear KGM_row Seevers_row SDR_row SOE_row
end




%% Plot the data
% figure(1)
% subplot(221)
% hold on
% %histogram(vertcat(SOE_meanKFactor{:,1}),50)
% histogram(vertcat(SDR_meanKFactor{:,1}),50)
% %histogram(vertcat(Seevers_meanKFactor{:,1}),50)
% histogram(vertcat(KGM_meanKFactor{:,1}),50)

%Wisconsin Data

%names = [{'Adams G5'} {'Adams G6'} {'Plainfield Lake PN1'} {'Plainfield Lake PN2'}];
names = [{'SDR'} {'Seevers'} {'SOE'} {'KGM'}];


for jj = 1:1
    
    figure(jj)

    hold on

    minDepth = min(depths{jj});
    maxDepth = max(depths{jj});


    subplot(141)
    % SDR 
    %plotIndex = plotIndex + 1;

    hold on
    for kk = 1:2000
        scatter(SDR_Kestimates{kk,jj},z{jj},4,'Filled','MarkerFaceColor',[152 152 152]./255)
    end
    scatter(SDR_Kmedian{jj},z{jj},20,'Filled')  %,'MarkerFaceColor',[1 0.2 0.2])
    scatter(K{jj}, z{jj},30)



    set(gca,'XScale','log')
    set(gca,'YDir','reverse')
    xlim([10^-6, 10^-2])
    ylim([1, 18])

    xlabel('K (m/s)')
    box on
    grid on
    set(gca,'FontSize',10)

    title(names(1))

    a=gca;
    a.XRuler.TickLabelGapOffset = -2;
    a.YRuler.TickLabelGapOffset = 1;

    set(gcf, 'Renderer', 'painters')
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Seevers
    subplot(142)
    hold on
    for kk = 1:2000
        scatter(Seevers_Kestimates{kk,jj},z{jj},4,'Filled','MarkerFaceColor',[152 152 152]./255)
    end
    scatter(Seevers_Kmedian{jj},z{jj},20,'Filled')  %,'MarkerFaceColor',[1 0.2 0.2])
    scatter(K{jj}, z{jj},30)



    set(gca,'XScale','log')
    set(gca,'YDir','reverse')
    xlim([10^-6, 10^-2])
    ylim([1, 18])

    xlabel('K (m/s)')
    box on
    grid on
    set(gca,'FontSize',10)

    title(names(2))

    a=gca;
    a.XRuler.TickLabelGapOffset = -2;
    a.YRuler.TickLabelGapOffset = 1;

    set(gcf, 'Renderer', 'painters')
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SOE
    subplot(143)
    hold on
    for kk = 1:1000
        scatter(SOE_Kestimates{kk,jj},z{jj},4,'Filled','MarkerFaceColor',[152 152 152]./255)
    end
    scatter(SOE_Kmedian{jj},z{jj},20,'Filled')  %,'MarkerFaceColor',[1 0.2 0.2])
    scatter(K{jj}, z{jj},30)



    set(gca,'XScale','log')
    set(gca,'YDir','reverse')
    xlim([10^-6, 10^-2])
    ylim([1, 18])

    xlabel('K (m/s)')
    box on
    grid on
    set(gca,'FontSize',10)

    title(names(3))

    a=gca;
    a.XRuler.TickLabelGapOffset = -2;
    a.YRuler.TickLabelGapOffset = 1;

    set(gcf, 'Renderer', 'painters')
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % KGM
    subplot(144)
    hold on
    for kk = 1:1000
        scatter(KGM_Kestimates{kk,jj},z{jj},4,'Filled','MarkerFaceColor',[152 152 152]./255)
    end
    scatter(KGM_Kmedian{jj},z{jj},20,'Filled')  %,'MarkerFaceColor',[1 0.2 0.2])
    scatter(K{jj}, z{jj},30)



    set(gca,'XScale','log')
    set(gca,'YDir','reverse')
    xlim([10^-6, 10^-2])
    ylim([1, 18])

    xlabel('K (m/s)')
    box on
    grid on
    set(gca,'FontSize',10)

    title(names(4))

    a=gca;
    a.XRuler.TickLabelGapOffset = -2;
    a.YRuler.TickLabelGapOffset = 1;

    set(gcf, 'Renderer', 'painters')


end

