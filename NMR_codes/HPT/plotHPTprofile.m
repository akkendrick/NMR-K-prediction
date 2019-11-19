% Plot HPT data profiles
% 
% load 'HPT_WellPN1.mat'
%tableName = PlainfieldHPT106052018NW;
% 
% load 'HPT_WellPN2.mat'
% tableName = PlainfieldHPT206052018SE;

% load 'HPT_WellG5'
% tableName = Adams05312018HPT1;

load 'HPT_WellG6'
tableName = Adams05312018HPT2B;
%%

depthft = tableName{:,1};
depthm = depthft * 0.3048;

EC = tableName{:,2};

figure(1)
hold on

plot(EC,depthm)

xlabel('EC (mS/m)')
ylabel('Depth (m)')
xlim([0 50])

grid on
box on

set(gca,'Ydir','reverse')