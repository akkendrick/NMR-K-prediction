% prospaT2.m
% Plot inverted CPMG data from Prospa in Matlab

% Alex Kendrick 
% 2/27/18

% Set inversion parameters
basicPath = strcat('/Users/akendrick/','Google Drive File Stream','/My Drive','/Stanford/Research/Research_Data/')
%basicPath = '/Users/akendrick/Google_Drive/Stanford/Research/Research_Data/';
%basicPath = 'G:\Google Drive\Stanford\Research\';
%basicPath = '/Users/akendrick/Desktop/';

% filePath = 'Research_Data/StCloud_multDelta/';
% experiment = 'DT2_wFlow_30000';
% dataSet = 'DT2_wFlow_30000/1Smooth/';
% filename = 'inversion.2d';

% filePath = 'Oct2018_glassBead/drained/';
% experiment = 'T2CPMG_glassBead_drained/';
% dataSet = 'inversion/';
% filename = 'spectrum.1d';
% 
% filePath = 'Sept2018_ZEOX/sample3/';
% experiment = 'T2CPMG_ZEOX_200/';
% dataSet = 'inversion/';
% filename = 'spectrum.1d';

filePath = 'Wisconsin_wellBulk/';
experiment = 'Plainfield_well2/T2CPMG/';
dataSet = "";
filename = 'spectrum.1d';



[x,y] = LoadProspaData(strcat(basicPath,filePath,experiment,dataSet,filename));
T2Params = LoadProspaParameters(strcat(basicPath,filePath,experiment,dataSet,'T2Analysis.par'));
timeVec = logspace(log10(T2Params.x_minimum), log10(T2Params.x_maximum), T2Params.x_steps);
%plot(timeVec, invertedData)

%yNorm = y ./ max(y);


%xcutoff = [33 33];
%ycutoff = [0 2];

semilogx(x,y, 'LineWidth',2)
xlabel('Relaxation Time T_2 (ms)')
ylabel('Relative Amplitude')
grid on 
box on

hold on

%semilogx(xcutoff, ycutoff)











