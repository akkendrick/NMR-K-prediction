% Analyze SDR b statistics (only designed for all data points currently)

sites = {'wisc_all','maurer_all'}
legendNames = {'Wisc','Maurer and Knight'}

% Load wisc and M&K data
for kk = 1:length(sites)
    %close all
    baseName = sites{kk}
    baseDir = '/Volumes/GoogleDrive/My Drive/USGS Project/USGS Data/';
    baseDir = '/Volumes/GoogleDrive/My Drive/USGS Project/NMR-K-prediction/Data/Aggregated_Data/';
    
    nmrName = baseName;
    
        [d, K, T2ML, phi, z, SumEch, logK, logT2ML, logPhi, SumEch_3s, SumEch_twm, ...
    SumEch_twm_3s] = loadnmrdata2(nmrName);

    SDR_K = @(b,m,n,phi,T2ML) (b.*(phi.^m).*(T2ML).^n);
    SDR_b = @(K,m,n,phi,T2ML) K./((phi.^m).*(T2ML).^n);
    
    SOE_K = @(C,SOE,n) C.*(SOE.^n);
    SOE_C = @(K,SOE,n) K./(SOE.^n);
    
    m = 0;
    n = 2;
    
    bProfile{kk} = SDR_b(K,m,n,phi,T2ML);
    cProfile{kk} = SOE_C(K, SumEch, n);
   
    KProfile{kk} = K;
    T2MLProfile{kk} = T2ML;
        
end
%%
% Load Neb data

load neb_data_3stk.mat

loggingDepth = neb_T2ML_stack3{:,1};

% Convert to seconds
T2ML = neb_T2ML_stack3{:,2}*10^-3;

WBF_topDepth = neb_WBF_K{:,1};
WBF_bottomDepth = neb_WBF_K{:,2};
WBF_K = neb_WBF_K{:,4}; % in ft/day

% Convert K data to m/s
WBF_K = WBF_K * (3.528*10^-6);

SDR_K = @(b,m,n,phi,T2ML) (b.*(phi.^m).*(T2ML).^n);
SDR_b = @(K,m,n,phi,T2ML) K./((phi.^m).*(T2ML).^n);

m = 0;
n = 2;

bstrpB = 0.0042;
color = [153 204 204]/255;
colorStar = [229 227 186]/255;

% Upscale NMR data to match WBF K intervals
for kk = 1:length(WBF_topDepth)
    currentDepths = loggingDepth(loggingDepth > WBF_topDepth(kk) &...
        loggingDepth < WBF_bottomDepth(kk));
    currentT2ML = T2ML(loggingDepth > WBF_topDepth(kk) &...
        loggingDepth < WBF_bottomDepth(kk));
    currentWBF_K = ones(length(currentDepths),1).*WBF_K(kk);
    
    currentSDRb = SDR_b(currentWBF_K,m,n,0.4,currentT2ML);
    
    %currentSDRK = SDR_K(
    
    % Follow upscaling equation from Dlubac et al. 2013
    p = length(currentDepths);
    
    intervalT2ML = 1/p*sum(currentT2ML);
    intervalb = 1/p*sum(currentSDRb);
    
    bNebProfile(kk) = intervalb;    
    T2MLNebProfile(kk) = intervalT2ML;
end

%%
KProfile{3} = WBF_K;
T2MLProfile{3} = T2MLNebProfile';
bProfile{3} = bNebProfile';


allK = vertcat(KProfile{:});
allb = vertcat(bProfile{:});

figure(2)
hold on
scatter(allK, allb,15,'filled')

grid on
box on

set(gca,'XScale','log')
set(gca,'YScale','log')

xlabel('K (m/s)')
ylabel('Corrected b')

%%
logAllK = log10(allK);

Nedges = 20;
[binsK, edgesK] = discretize(logAllK,Nedges); 

figure()
h = histogram(logAllK,edgesK)
%binnedSDRb = discretize(allb,Nedges

for jj = 1:length(edgesK)
    binnedK{jj} = allK(binsK == jj);
    binnedlogK{jj} = logAllK(binsK == jj);
    binnedb{jj} = allb(binsK == jj);
    
    if ~isempty(binnedb{jj})
        meanbinb(jj) = mean(binnedb{jj});
        maxbinb(jj) = max(binnedb{jj});
        minbinb(jj) = min(binnedb{jj});
        medianbinb(jj) = median(binnedb{jj});

        meanbinK(jj) = mean(binnedK{jj});

        meanbinLogK(jj) = mean(binnedlogK{jj});
    else
        meanbinb(jj) = NaN;
        maxbinb(jj) = NaN;
        minbinb(jj) =NaN;
        medianbinb(jj) = NaN;

        meanbinK(jj) = NaN;

        meanbinLogK(jj) = NaN;
    end
end

blogDiff = log10(maxbinb) - log10(minbinb)
bDiff = 10.^blogDiff

figure()
hold on

scatter(10.^meanbinLogK, meanbinb, 15,'filled')
scatter(10.^meanbinLogK, maxbinb, 15,'filled')
scatter(10.^meanbinLogK, minbinb, 15,'filled')

set(gca,'XScale','log')
set(gca,'YScale','log')

grid on
box on


    