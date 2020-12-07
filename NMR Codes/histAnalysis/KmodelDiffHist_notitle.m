function KmodelDiffHist_notitle(histDataSDR,histDataSOE, histDataSeevers, histDataKGM)

% Plot histogram of Kdiff Factor
%edges = [1 1.5 2 2.5 3 4 5 10 20 40 80];

edgesLog = [-2 -1.5 -1 -0.8 -0.6 -0.4 -0.2 0 0.2 0.4 0.6 0.8 1 1.5 2];
edges = 10.^edgesLog;


subplot(2,2,1)
hold on
grid on
box on
grid minor

%histogram(SDR_diffFactor_all, edges)
histogram(histDataSDR,edgesLog)
%histogram(neghistSDR,edgesLog)


ylim([0 50])
%ylim([0 30])

xticklabels({'-100','-10','0','10','100'})
xlabel('K Difference Factor')
ylabel('Counts')
set(gca,'FontSize',14)

subplot(2,2,4)
hold on
grid on
box on
grid minor

%histogram(Seevers_diffFactor_all, edges)
histogram(histDataSeevers,edgesLog)
%histogram(neghistSeevers,edgesLog)

ylim([0 50])
%ylim([0 30])

xticklabels({'-100','-10','0','10','100'})
xlabel('K Difference Factor')
set(gca,'FontSize',14)


subplot(2,2,3)

hold on
grid on
box on
grid minor

%histogram(SOE_diffFactor_all, edges)
histogram(histDataSOE,edgesLog)
%histogram(neghistSOE,edgesLog)

ylim([0 50])
%ylim([0 30])

xticklabels({'-100','-10','0','10','100'})
xlabel('K Difference Factor')
ylabel('Counts')
set(gca,'FontSize',14)


subplot(2,2,2)
hold on
grid on
box on
grid minor

%histogram(KGM_diffFactor_all, edges)
histogram(histDataKGM,edgesLog,'EdgeColor',[0 0.4470 0.7410],'FaceColor',[0 0.4470 0.7410],'FaceAlpha',0.6)
%histogram(neghistKGM,edgesLog,'EdgeColor',[0 0.4470 0.7410],'FaceColor',[0 0.4470 0.7410],'FaceAlpha',0.6)

ylim([0 50])
%ylim([0 25])

xticklabels({'-100','-10','0','10','100'})
xlabel('K Difference Factor')
set(gca,'FontSize',14)


end