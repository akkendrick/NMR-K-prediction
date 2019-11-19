% Load and plot nebraska T2/K data
clear
load neb_data_3stk.mat

loggingDepth = neb_T2ML_stack3{:,1};

% Convert to seconds
T2ML = neb_T2ML_stack3{:,2}*10^-3;

WBF_topDepth = neb_WBF_K{:,1};
WBF_bottomDepth = neb_WBF_K{:,2};
WBF_K = neb_WBF_K{:,4}; % in ft/day

% Convert K data to m/s
WBF_K = WBF_K * (3.528*10^-6);

SDR_K = @(b,m,n,phi,T2ML) (b.*(phi.^m).*(T2ML).^n);
SDR_b = @(K,m,n,phi,T2ML) K./((phi.^m).*(T2ML).^n);

m = 0;
n = 1;

bstrpB = 0.0042;
color = [153 204 204]/255;
colorStar = [229 227 186]/255;

% Upscale NMR data to match WBF K intervals
for kk = 1:length(WBF_topDepth)
    currentDepths = loggingDepth(loggingDepth > WBF_topDepth(kk) &...
        loggingDepth < WBF_bottomDepth(kk));
    currentT2ML = T2ML(loggingDepth > WBF_topDepth(kk) &...
        loggingDepth < WBF_bottomDepth(kk));
    currentWBF_K = ones(length(currentDepths),1).*WBF_K(kk);
    
    currentSDRb = SDR_b(currentWBF_K,m,n,0.4,currentT2ML);
    
    %currentSDRK = SDR_K(
    
    % Follow upscaling equation from Dlubac et al. 2013
    p = length(currentDepths);
    
    intervalT2ML = 1/p*sum(currentT2ML);
    intervalb = 1/p*sum(currentSDRb);
    
    bProfile(kk) = intervalb;    
    T2MLProfile(kk) = intervalT2ML;
end

figure(6)
hold on

scatter(log10(T2MLProfile),log10(WBF_K),30,'Filled')


grid on 
box on
% set(gca,'XScale','log')
% set(gca,'YScale','log')
% xlabel('T2ML (s)')
% ylabel('WBF K (m/s)')


KLog = log10(WBF_K);
bLog = log10(bProfile)';

% Compute linear fit
X = [ones(length(KLog),1) KLog];
fitParams = X\bLog;

slope = fitParams(2)
yInt = fitParams(1)


figure(2)
hold on

logMeanK = 10.^mean(log10(WBF_K));

scatter(WBF_K,bProfile,25,'Filled','MarkerEdgeColor',color, 'MarkerFaceColor',color)
plot(logMeanK,bstrpB,'^','MarkerSize',15,'LineWidth',2,'MarkerEdgeColor','k')

%plot(logMeanK,bstrpB,'p','MarkerSize',20,'MarkerFaceColor',colorStar,'MarkerEdgeColor','k')

grid on
box on

set(gca,'XScale','log')
set(gca,'YScale','log')

xlabel('K (m/s)')
ylabel('Corrected b')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(3)
hold on

logMeanT2ML = 10.^mean(log10(T2MLProfile));

scatter(T2MLProfile,bProfile,15,'Filled','MarkerEdgeColor',color, 'MarkerFaceColor',color)
plot(logMeanT2ML,bstrpB,'p','MarkerSize',20,'MarkerFaceColor',colorStar,'MarkerEdgeColor','k')

