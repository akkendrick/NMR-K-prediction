% Estimate error factors

load bestKModels.mat

kModel = dlubacModel;

for kk = 1:length(dlubacModel)
    errorEst(kk,:) = computeError(
