% Investigate CPMG Data

site = {'Site1-WellG6'}

[T2dist,T2logbins,SEdecayTime,SEdecayUniform,SEdecay,oneDVectors,...
    oneDVectorsUniform, nmrName] = loadAllRawNMRdata(site);

depths = SEdecay(:,1);
SEdecay = SEdecay(:,2:end);
SEdecayUniform = SEdecayUniform(:,2:end);

slice = 30;

figure()
hold on
plot(SEdecayTime, SEdecay(slice,:))
%plot(SEdecayTime, SEdecayUniform(slice,:))