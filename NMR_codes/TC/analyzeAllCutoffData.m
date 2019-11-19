% Analyze cutoff data 
% Script to look at optimal cutoff data and figure out what is going on
clear

load optimalCutoffTable_n2_m0.mat
%goodSites = [2,5,3,6,7,8];
goodSites = [1,2,3,4];
%goodSites = [1,2,3,4];
figure(1)

cutoff = cutoff*10^3;

for kk = 1:length(goodSites)

    index = goodSites(kk);
    
    smoothWindow = 20;
    smoothCurve = smooth(totalErrorMatrix(index,:),smoothWindow);
    [minVal, minIndex] = min(smoothCurve);
    
    optimalCutoff(kk) = cutoff(minIndex);
    
    subplot(2,2,kk)
    %subplot(3,2,kk)
    hold on
    grid on
    box on
    
    plot(cutoff,smoothCurve,'HandleVisibility','off','LineWidth',2)
    scatter(cutoff,totalErrorMatrix(index,:),5,'Filled')

    plot(cutoff(minIndex),smoothCurve(minIndex),'k*','MarkerSize',10,'HandleVisibility','off')
    
    ylim([0.0,0.1])
   % xlim([0.02,0.5])
   xlim([10,500]) 
   title(siteList(index))
   xlabel('Cutoff (ms)')
   ylabel('MAE relative to K_{DPP}')
    
    
end

%legend('n1 m1','n1 m2','n1 mvar','n2 m1','n2 m2','n2 mvar')
%legend('n1 m1','n2 m1')

%%
load optimalCutoffTableWash_n2_m0.mat
%goodSites = [2,5,3,6,7,8];
goodSites = [1,2];
%goodSites = [1,2,3,4];
figure(1)

cutoff = cutoff*10^3;

for kk = 1:length(goodSites)

    index = goodSites(kk);
    
    smoothWindow = 20;
    smoothCurve = smooth(totalErrorMatrix(index,:),smoothWindow);
    [minVal, minIndex] = min(smoothCurve);
    optimalCutoff(kk) = cutoff(minIndex);

    subplot(2,2,kk)
    %subplot(3,2,kk)
    hold on
    grid on
    box on
    
    plot(cutoff,smoothCurve,'HandleVisibility','off','LineWidth',2)
    scatter(cutoff,totalErrorMatrix(index,:),5,'Filled')

    plot(cutoff(minIndex),smoothCurve(minIndex),'k*','MarkerSize',10,'HandleVisibility','off')
    
    ylim([0.0,0.1])
   % xlim([0.02,0.5])
   xlim([10,500]) 
   title(siteList(index))
   xlabel('Cutoff (ms)')
   ylabel('MAE relative to K_{DPP}')
    
    
end

%%
load optimalCutoffTableKansas_n2_m0.mat
%goodSites = [2,5,3,6,7,8];
goodSites = [1,2,3,4,5,6,7,8];
%goodSites = [1,2,3,4];
figure(1)

cutoff = cutoff*10^3;

for kk = 1:length(goodSites)

    index = goodSites(kk);
    
    smoothWindow = 20;
    smoothCurve = smooth(totalErrorMatrix(index,:),smoothWindow);
    [minVal, minIndex] = min(smoothCurve);
    
    optimalCutoff(kk) = cutoff(minIndex);

    
    subplot(4,2,kk)
    %subplot(3,2,kk)
    hold on
    grid on
    box on
    
    plot(cutoff,smoothCurve,'HandleVisibility','off','LineWidth',2)
    scatter(cutoff,totalErrorMatrix(index,:),5,'Filled')

    plot(cutoff(minIndex),smoothCurve(minIndex),'k*','MarkerSize',10,'HandleVisibility','off')
    
    ylim([0.0,0.1])
   % xlim([0.02,0.5])
   xlim([10,500]) 
   title(siteList(index))
   xlabel('Cutoff (ms)')
   ylabel('MAE relative to K_{DPP}')
    
    
end