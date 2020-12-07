clear
load ostermanClay.mat

% Take homogeneous data for now
clayPercent = OstermanclaySDRbK{1:12,1};
K = OstermanclaySDRbK{1:12,2};
T2ML = OstermanclaySDRbK{1:12,3};
phi = OstermanclaySDRbK{1:12,4};

SDR_b = @(K,m,n,phi,T2ML) K./((phi.^m).*(T2ML).^n);

bProfile = SDR_b(K,0,2,phi,T2ML);

logMeanK = 10.^mean(log10(K));


figure(2)
hold on

scatter(K, bProfile,25,'Filled','r') 
%plot(logMeanK, 4.71*10^-2,'p','MarkerSize',25,'MarkerFaceColor','r','MarkerEdgeColor','k')


grid on
box on

set(gca,'XScale','log')
set(gca,'YScale','log')

xlabel('K (m/s)')
ylabel('Corrected b')

% Look at clustered Data
clayPercentClust = OstermanclaySDRbK{14:22,1};
KClust = OstermanclaySDRbK{14:22,2};
T2MLClust = OstermanclaySDRbK{14:22,3};
phiClust = OstermanclaySDRbK{14:22,4};

bProfileClust = SDR_b(KClust,0,2,phiClust,T2MLClust);

figure(2)
scatter(KClust, bProfileClust)