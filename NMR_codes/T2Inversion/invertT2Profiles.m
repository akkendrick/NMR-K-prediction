% Invert Dart T2 profiles
close all
clear


site = {'Site2-WellPN2'}

[T2dist,T2logbins,SEdecayTime,SEdecayUniform,SEdecay,oneDVectors,...
    oneDVectorsUniform, nmrName] = loadAllRawNMRdata(site);

depths = T2dist(:,1);
SEdecay = SEdecay(:,2:end);
SEdecayUniform = SEdecayUniform(:,2:end);
noise = oneDVectors(:,14);
porosity = oneDVectors(:,2);
T2dist = T2dist(:,2:end);

T2linbins = 10.^T2logbins;

newT2linbins = logspace(-3.1037, 0.6990, 200);
newT2logbins = log10(newT2linbins);
%% 
% slice = 30;
% smoothing = 40;
% 
% I = unifrnd(-noise(slice),noise(slice),size(SEdecay,2),1)';
% [m,dsyn] = whittallmodt2fit(SEdecay(slice,:),SEdecayTime,I,newT2linbins,smoothing);
% 
% 
% figure(1)
% hold on
% plot(SEdecayTime, SEdecay(slice,:))
% plot(SEdecayTime, I)
% plot(dsyn(:,1), dsyn(:,2))
% plot(dsyn(:,1), dsyn(:,3))
% 
% figure(2)
% hold on
% plot(T2linbins, T2dist(slice,:)./max(T2dist(slice,:)))
% plot(newT2linbins, m(1:end-1)./max(m(1:end-1)))
% 
% set(gca,'XScale','log')
%%
% Invert the entire profile
smoothing = 0.1;
nBins = length(newT2linbins);
nDepths = length(depths)
newT2data = zeros(nDepths,nBins);
for jj = 1:length(depths)
    slice = jj;
    I = unifrnd(-noise(slice),noise(slice),size(SEdecay,2),1)';
    [m,dsyn] = whittallmodt2fit(SEdecay(slice,:),SEdecayTime,I,newT2linbins,smoothing);

    scaledT2data = normalize(m(1:nBins),'norm',1);
    %scaledT2data = m(1:nBins);
    newT2data(jj,1:nBins) = scaledT2data';
    
    simulatedDecayTime(jj,:) = dsyn(:,1);
    simulatedImagDecay(jj,:) = dsyn(:,3);
    
    T2lm_new(jj,1) = 10.^(sum(newT2logbins.*newT2data(jj,1:nBins))./(sum(newT2data(jj,1:nBins))));
    T2lm(jj,1) = 10.^(sum(T2logbins.*T2dist(jj,:))./porosity(jj));

    
end

%%
slice = 30;

figure(1)
hold on
T2sliceNorm = normalize(T2dist(slice,:),'norm',1);
plot(T2linbins, T2sliceNorm)
plot(newT2linbins, newT2data(slice,:))
plot(T2lm_new(slice),0.05,'r*')
plot(T2lm(slice),0.05,'g*')
set(gca,'XScale','log')

figure(2)
hold on
plot(T2lm_new, 'r')
plot(T2lm, 'g')


%%
newT2Dist = [depths newT2data];

if strcmp(site,'Site1-WellG5')
    name = 'G5_W1_tr5_20x_16p5_up_F1n2_wRIN_wRFI_Reg50_Va1';
elseif strcmp(site,'Site1-WellG6')
    name = 'G6_W2_tr5_20x_16p75_up_F_wRIN_wRFI_reg50_Va1';
elseif strcmp(site,'Site2-WellPN1')
    name = 'Pl_W1_Tr5_20x_MPp75aLS_F1n2_wRIN_wRFI_Reg50_Va1';
elseif strcmp(site,'Site2-WellPN2')
    name = 'W2_Tr5_20x_MPp75aLS_Reg50_wRIN_wRFI_Va1';
end
outFileName = [name '_' num2str(smoothing) '_newInvert.mat'];
save(outFileName,'newT2Dist','newT2linbins')


