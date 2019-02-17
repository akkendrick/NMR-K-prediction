% Run Timur-Coates Model Estimates


% Range over pairs of m and n values
close all
clear

load enso 

siteList = [{'Site1-WellG6above'} {'Site1-WellG6below'} {'Site1-WellG5above'}...
    {'Site1-WellG5below'} {'Site2-WellPN1'} {'Site2-WellPN2'} {'Site1-WellG6'} {'Site1-WellG5'}];

% Matrix Structure... Repeat for each site
%       Bootstrap           
% c
% m
% n

%cutoff = [20 33 40 50 60 80 100 120 140 160 200 240 280 320 380 440 500 560 620]*10^-3;
cutoff = (30:10:800)*10^-3;
%cutoff = 510*10^-3;
m = 2;
n = 1;

figureson = 0;
wDirect = 1;
currentRow = -3;
%%
parfor j = 1:length(cutoff)
    cutoff(j)
    currentFitMatrix = [];
    currentErrorMatrix = [];

    for i = 1:length(siteList)
        baseName = siteList{i}

        [K,z,T2dist,T2logbins,kTC_best,bestFitMatrix,totalErrorEstimate] = computeTCperm(baseName,n,m,cutoff(j),figureson);

        currentFitMatrix = [currentFitMatrix bestFitMatrix];
        currentErrorMatrix = [currentErrorMatrix totalErrorEstimate];
    end

    totalErrorMatrix(j,:) = currentErrorMatrix;
end
        
save('optimalCutoffTable.mat','totalErrorMatrix','cutoff','siteList','n','m')

%%
load('optimalCutoffTable.mat')
cutoff = (20:10:800)*10^-3;

g6aSmooth = smooth(totalErrorMatrix(:,1),10);
g5aSmooth = smooth(totalErrorMatrix(:,3),10);
g6bSmooth = smooth(totalErrorMatrix(:,2),10);
g5bSmooth = smooth(totalErrorMatrix(:,4),10);
pn1Smooth = smooth(totalErrorMatrix(:,5),10);
pn2Smooth = smooth(totalErrorMatrix(:,6),10); 
g6Smooth = smooth(totalErrorMatrix(:,7),10);
g5Smooth = smooth(totalErrorMatrix(:,8),10);

smoothErrorMatrix = [g6aSmooth g6bSmooth g5aSmooth g5bSmooth pn1Smooth pn2Smooth g6Smooth g5Smooth];

[minVals, index] = min(smoothErrorMatrix);
%g6aFit = fit(cutoff(15:end)',totalErrorMatrix(15:end,1),'smoothingspline','SmoothingParam',0.99);
cutoff_ms = cutoff*10^3;


figure(1)
subplot(2,2,1)
hold on
plot(cutoff_ms, [totalErrorMatrix(:,1) totalErrorMatrix(:,3)], '.')
%plot(smoothCutoff, g6aFit(smoothCutoff))
plot(cutoff_ms, g6aSmooth,'LineWidth',2)
plot(cutoff_ms, g5aSmooth,'LineWidth',2)
plot([cutoff_ms(index(1)) cutoff_ms(index(3))],[smoothErrorMatrix(index(1),1) smoothErrorMatrix(index(3),3)],'k*','MarkerSize',10)

legend([siteList{1,1}; siteList{1,3}])
grid on
box on
xlabel('Cutoff (ms)')
ylabel('Mean Average Error (m/day)')
ylim([0.0002,0.001])

subplot(2,2,2)
hold on
plot(cutoff_ms, [totalErrorMatrix(:,2) totalErrorMatrix(:,4)], '.')
plot(cutoff_ms, g6bSmooth,'LineWidth',2)
plot(cutoff_ms, g5bSmooth,'LineWidth',2)
plot([cutoff_ms(index(2)) cutoff_ms(index(4))],[smoothErrorMatrix(index(2),2) smoothErrorMatrix(index(4),4)],'k*','MarkerSize',10)

legend([siteList{1,2}; siteList{1,4}])
grid on
box on
xlabel('Cutoff (ms)')
ylabel('Mean Average Error (m/day)')
ylim([0.001,0.005])

subplot(2,2,3)
hold on
plot(cutoff_ms, [totalErrorMatrix(:,5) totalErrorMatrix(:,6)], '.')
plot(cutoff_ms, pn1Smooth,'LineWidth',2)
plot(cutoff_ms, pn2Smooth,'LineWidth',2)
plot([cutoff_ms(index(5)) cutoff_ms(index(6))],[smoothErrorMatrix(index(5),5) smoothErrorMatrix(index(6),6)],'k*','MarkerSize',10)

legend([siteList{1,5}; siteList{1,6}])
grid on
box on
xlabel('Cutoff (ms)')
ylabel('Mean Average Error (m/day)')
ylim([0.0015,0.005])

subplot(2,2,4)
hold on
plot(cutoff_ms, [totalErrorMatrix(:,7) totalErrorMatrix(:,8)], '.')
plot(cutoff_ms, g6Smooth,'LineWidth',2)
plot(cutoff_ms, g5Smooth,'LineWidth',2)
plot([cutoff_ms(index(7)) cutoff_ms(index(8))],[smoothErrorMatrix(index(7),7) smoothErrorMatrix(index(8),8)],'k*','MarkerSize',10)


legend([siteList{1,7}; siteList{1,8}])
grid on
box on
xlabel('Cutoff (ms)')
ylabel('Mean Average Error (m/day)')
ylim([0.0015,0.005])



