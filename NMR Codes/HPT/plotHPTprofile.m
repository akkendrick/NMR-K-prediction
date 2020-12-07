% Plot HPT data profiles
% % 
% load 'HPT_WellPN1.mat'
%tableName = PlainfieldHPT106052018NW;
% 
% load 'HPT_WellPN2.mat'
% tableName = PlainfieldHPT206052018SE;

% load 'HPT_WellG5'
% tableName = Adams05312018HPT1;
% % 
load 'HPT_WellG6'
tableName = Adams05312018HPT2B;
%%

depthft = tableName{:,1};
%depthm = depthft * 0.3048;
depthm = tableName{:,13};

EC = tableName{:,2};

figure(5)
hold on

plot(EC,depthm)

xlabel('EC (mS/m)')
ylabel('Depth (m)')
xlim([0 50])

grid on
box on

set(gca,'Ydir','reverse')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%

HPT_pressAvgPSI = tableName{:,5};

figure(2)
hold on

plot(HPT_pressAvgPSI,depthm,'LineWidth',2)

xlabel('HPT Avg Pressure (psi)')
ylabel('Depth (m)')
xlim([0 60])
ylim([0 11])

legend('G5','G6')

grid on
box on

set(gca,'Ydir','reverse')
set(gca,'FontSize',12)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

HPT_flowAvgmLmin = tableName{:,8};

figure(3)
hold on

plot(HPT_flowAvgmLmin,depthm)

xlabel('HPT Avg Flow (ml/min)')
ylabel('Depth (m)')
%xlim([0 50])

grid on
box on

set(gca,'Ydir','reverse')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%

HPT_linePress_Avg = tableName{:,11};

figure(4)
hold on

plot(HPT_linePress_Avg,depthm)

xlabel('HPT Avg Line Pressure (psi)')
ylabel('Depth (m)')
%xlim([0 50])

grid on
box on

set(gca,'Ydir','reverse')





