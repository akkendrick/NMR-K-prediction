function KmodelDiffHist_TC(histDataSDR,histData_TCref, histData_TCideal)

% Plot histogram of Kdiff Factor
%edges = [1 1.5 2 2.5 3 4 5 10 20 40 80];

edgesLog = [-2 -1.5 -1 -0.8 -0.6 -0.4 -0.2 0 0.2 0.4 0.6 0.8 1 1.5 2];
edges = 10.^edgesLog;


subplot(3,1,1)

title('SDR Model')
hold on
grid on
box on
grid minor

%histogram(SDR_diffFactor_all, edges)
histogram(histDataSDR,edgesLog)
%histogram(neghistSDR,edgesLog)


ylim([0 50])
%ylim([0 30])

xticklabels({'-100','-32','-10','-0.32','0','0.32','10','32','100'})
xlabel('K Difference Factor')
ylabel('Counts')
set(gca,'FontSize',14)

subplot(3,1,2)
title('Median min T_2 from bimodal plots cutoff')
hold on
grid on
box on
grid minor

%histogram(Seevers_diffFactor_all, edges)
histogram(histData_TCideal,edgesLog)
%histogram(neghistSeevers,edgesLog)

ylim([0 50])
%ylim([0 30])

xticklabels({'-100','-32','-10','-0.32','0','0.32','10','32','100'})
%xticklabels({'-100','-10','0','10','100'})
xlabel('K Difference Factor')
set(gca,'FontSize',14)


subplot(3,1,3)

title('33 ms cutoff')
hold on
grid on
box on
grid minor

%histogram(SOE_diffFactor_all, edges)
histogram(histData_TCref,edgesLog)
%histogram(neghistSOE,edgesLog)

ylim([0 50])
%ylim([0 30])

xlabel('K Difference Factor')
ylabel('Counts')
set(gca,'FontSize',14)

xticklabels({'-100','-32','-10','-0.32','0','0.32','10','32','100'})
%xticklabels({'-100','-10','0','10','100'})
xlabel('K Difference Factor')
set(gca,'FontSize',14)


end