% Run Timur-Coates Model Estimates


% Range over pairs of m and n values
close all
clear

load enso 

siteList = [{'Site1-WellG6'} {'Site1-WellG6above'} {'Site1-WellG6below'} {'Site1-WellG5'} {'Site1-WellG5above'}...
    {'Site1-WellG5below'} {'Site2-WellPN1'} {'Site2-WellPN2'}  ];

% Matrix Structure... Repeat for each site
%       Bootstrap           
% c
% m
% n

%cutoff = [20 33 40 50 60 80 100 120 140 160 200 240 280 320 380 440 500 560 620]*10^-3;
 cutoff_1 = (10:2:800)*10^-3;
 cutoff_2 = (800:20:2200)*10^-3;
 cutoff = [cutoff_1 cutoff_2];
%cutoff = [20 200];

%cutoff = 510*10^-3;
m = [1];
n = [1];

figureson = 0;
wDirect = 1;
currentRow = -3;
%%
tic
totalErrorMatrix = zeros(length(siteList),length(cutoff));
totalmMatrix = zeros(length(siteList),length(cutoff));
totalnMatrix = zeros(length(siteList),length(cutoff));
totalcMatrix = zeros(length(siteList),length(cutoff));

for j = 1:length(siteList)
    currentFitMatrix = [];
    currentErrorMatrix = [];
    baseName = siteList{j};

    mTemp = zeros(1,length(cutoff));
    nTemp = zeros(1,length(cutoff));
    cTemp = zeros(1,length(cutoff));
    errorTemp = zeros(1,length(cutoff));
    
    
    parfor i = 1:length(cutoff)

        [K,z,T2dist,T2logbins,kTC_best,bestFitMatrix,totalErrorEstimate] = computeTCperm(baseName,n,m,cutoff(i),figureson);

        mTemp(i) = bestFitMatrix(3);
        nTemp(i) = bestFitMatrix(2);
        cTemp(i) = bestFitMatrix(1);
        errorTemp(i) = totalErrorEstimate;
    end
    
    totalmMatrix(j,:) = mTemp;
    totalnMatrix(j,:) = nTemp;
    totalcMatrix(j,:) = cTemp;
    totalErrorMatrix(j,:) = errorTemp;
        
end
toc

save('optimalCutoffTable_n1_m1.mat','totalmMatrix','totalnMatrix','totalcMatrix','totalErrorMatrix','cutoff','siteList','n','m')

%%
load('optimalCutoffTable_n1_m1.mat')
%cutoff = (20:2:800)*10^-3;

smoothWindow = 40;

siteList = [{'Site1-WellG6'} {'Site1-WellG6above'} {'Site1-WellG6below'} {'Site1-WellG5'} {'Site1-WellG5above'}...
    {'Site1-WellG5below'} {'Site2-WellPN1'} {'Site2-WellPN2'}  ];

g6Smooth = smooth(totalErrorMatrix(1,:),smoothWindow);
g5Smooth = smooth(totalErrorMatrix(4,:),smoothWindow);

g6aSmooth = smooth(totalErrorMatrix(2,:),smoothWindow);
g5aSmooth = smooth(totalErrorMatrix(5,:),smoothWindow);

g6bSmooth = smooth(totalErrorMatrix(3,:),smoothWindow);
g5bSmooth = smooth(totalErrorMatrix(6,:),smoothWindow);

pn1Smooth = smooth(totalErrorMatrix(7,:),smoothWindow);
pn2Smooth = smooth(totalErrorMatrix(8,:),smoothWindow); 

smoothErrorMatrix = [g6aSmooth g6bSmooth g5aSmooth g5bSmooth pn1Smooth pn2Smooth g6Smooth g5Smooth]';

[minVals, index] = min(smoothErrorMatrix,[],2);
%g6aFit = fit(cutoff(15:end)',totalErrorMatrix(15:end,1),'smoothingspline','SmoothingParam',0.99);
cutoff_ms = cutoff*10^3;


figure(1)
subplot(2,2,1)
hold on
plot(cutoff_ms, [totalErrorMatrix(2,:); totalErrorMatrix(5,:)], '.')
%plot(smoothCutoff, g6aFit(smoothCutoff))
plot(cutoff_ms, g6aSmooth,'LineWidth',2)
plot(cutoff_ms, g5aSmooth,'LineWidth',2)
plot([cutoff_ms(index(1)) cutoff_ms(index(3))],[smoothErrorMatrix(1,index(1)); smoothErrorMatrix(3,index(3))],'k*','MarkerSize',10)

legend({siteList{1,2}; siteList{1,5}})
grid on
box on
xlabel('Cutoff (ms)')
ylabel('Mean Average Error (m/day)')
%ylim([0.0002,0.001])
ylim([0,1])

subplot(2,2,2)
hold on
plot(cutoff_ms, [totalErrorMatrix(3,:); totalErrorMatrix(6,:)], '.')
plot(cutoff_ms, g6bSmooth,'LineWidth',2)
plot(cutoff_ms, g5bSmooth,'LineWidth',2)
plot([cutoff_ms(index(2)) cutoff_ms(index(4))],[smoothErrorMatrix(2,index(2)); smoothErrorMatrix(4,index(4))],'k*','MarkerSize',10)

legend({siteList{1,3}; siteList{1,6}})
grid on
box on
xlabel('Cutoff (ms)')
ylabel('Mean Average Error (m/day)')
%ylim([0.001,0.005])
ylim([0,1])

subplot(2,2,3)
hold on
plot(cutoff_ms, [totalErrorMatrix(7,:); totalErrorMatrix(8,:)], '.')
plot(cutoff_ms, pn1Smooth,'LineWidth',2)
plot(cutoff_ms, pn2Smooth,'LineWidth',2)
plot([cutoff_ms(index(5)) cutoff_ms(index(6))],[smoothErrorMatrix(5,index(5)); smoothErrorMatrix(6,index(6))],'k*','MarkerSize',10)

legend({siteList{1,7}; siteList{1,8}})
grid on
box on
xlabel('Cutoff (ms)')
ylabel('Mean Average Error (m/day)')
%ylim([0.0015,0.005])
ylim([0,1])

subplot(2,2,4)
hold on
plot(cutoff_ms, [totalErrorMatrix(1,:); totalErrorMatrix(4,:)], '.')
plot(cutoff_ms, g6Smooth,'LineWidth',2)
plot(cutoff_ms, g5Smooth,'LineWidth',2)
plot([cutoff_ms(index(7)) cutoff_ms(index(8))],[smoothErrorMatrix(7,index(7)); smoothErrorMatrix(8,index(8))],'k*','MarkerSize',10)


legend({siteList{1,1}; siteList{1,4}})
grid on
box on
xlabel('Cutoff (ms)')
ylabel('Mean Average Error (m/day)')
%ylim([0.0015,0.005])
ylim([0,1])



