% Plot gamma/EMI logs
close all

baseDir = '/Volumes/GoogleDrive/My Drive/Stanford/USGS Project/Field Data/USGS Data/WI_gamma-EMI-bLS_csvfiles';

siteName = 'Site1-G6-Well2-gamma-EMI-bLS.csv'

in1 = [baseDir '/' siteName];

gammaEMIdata = csvread(in1,3,0);

depth = gammaEMIdata(:,1);
gamma = gammaEMIdata(:,2);
EMI = gammaEMIdata(:,3);

EMIdepth = depth(EMI > 0);
goodEMI = EMI(EMI > 0);

figure(1)
subplot(121)
hold on
z = ones(1,length(depth));

% plot3(gamma, depth, z)
% area(gamma)
plot(gamma, depth)

%view(0,90)
set(gca,'YDir','reverse')

grid on
box on

subplot(122)


plot(goodEMI, EMIdepth)
set(gca,'YDir','reverse')

grid on
box on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(2)
ax1 = axes('Position',[0.1 0.1 0.8 0.8]);
hold(ax1, 'on')

smoothWindow = 10;
smoothGamma = smooth(gamma, smoothWindow);

%plot(gamma, depth, 'parent', ax1)
plot(smoothGamma, depth,'parent',ax1,'LineWidth',2)

ylim([0,18])
set(ax1,'YDir','reverse')

ax1.XColor = 'k';
ax1.YColor = 'k';

xlabel('Counts/second')
ylabel('Depth (m)')

box on 
grid on

ax1_pos = ax1.Position; % position of first axes
ax2 = axes('Position',ax1_pos,...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','none','NextPlot','add');
ax2.XColor = 'r';

hold(ax2, 'on')


plot(10+goodEMI, EMIdepth, 'r', 'parent', ax2,'LineWidth',2)
set(gca,'YDir','reverse')
set(gca,'YTickLabel',[])
ylim([0,18])
xlim([0,90])

xlabel('10+mS/m')
