%Plot NMR noise for each well

sites = [{'Site1-WellG5'} {'Site1-WellG6'} {'Site2-WellPN1'} {'Site2-WellPN2'}];

for kk = 1:length(sites)
    [T2dist, T2logbins, SEdecayTime, SEdecayUniform, SEdecay,...
        oneDVectors, oneDVectorsUniform,nmrName] = loadAllRawNMRdata(sites{kk});
   
    [d, K, T2ML, phi, z, SumEch, logK, logT2ML, logPhi, SumEch_3s, SumEch_twm, ...
    SumEch_twm_3s] = loadnmrdata2(nmrName);
    
    Kall{kk} = K;
    Phiall{kk} = phi;
    T2MLall{kk} = T2ML;
    SumEchAll{kk} = SumEch;
    noise{kk} = oneDVectors(:,14);
   
end

load('G7_data.mat')
noise{5} = G7Tr530x1p0F2wRINnoRFIreg50Va1.Noise;

names = [{'Site1-WellG5'} {'Site1-WellG6'} {'Site2-WellPN1'} {'Site2-WellPN2'} {'Site1-WellG7'}];
edges = [2 2.5 3 3.5 4 4.5 5 5.5 6 6.5 7 8 9 10 12 ];


for kk = 1:length(sites) + 1
   subplot(3,2,kk)
   
   histogram(noise{kk},edges)
   grid on
   box on
   
   title(names{kk})
   xlabel('NMR Noise')
   ylabel('Number of Counts')
      
   
end