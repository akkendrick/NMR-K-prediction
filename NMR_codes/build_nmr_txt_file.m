function [z, T2ML, phi, Dk, soe, soe_3s, soe_twm, soe_twm_3s,k_SDR_VC, t, S, CapH20, FreeH20, k_TC] = build_nmr_txt_file(site,name,offset)
% convert raw data files to text files in format to be used by
% NMR_script_main

% read in data files
% in1 = [pwd '\bstrp_dat_VC_raw\' name '\' name '_SE_decay' '.txt']; 
% in2 = [pwd '\bstrp_dat_VC_raw\' name '\' name '.txt'];
% in3 = [pwd '\bstrp_dat_original\' name '.txt'];



baseDir = '/Volumes/GoogleDrive/My Drive/USGS Project/USGS Data/';

in1 = [baseDir site '/' name '/' name '_SE_decay' '.txt']; 
in2 = [baseDir site '/' name '/' name '_1Dvectors.txt'];
in3 = [baseDir site '/' name '/' strcat(site,'_DPP_filt.txt')];
in4 = [baseDir site '/' name '/' name '_SE_decay_time.txt'];

decaycurve = load(in1); 
dparam = dlmread(in2,'\t',1,0); 
DPPdat = load(in3); 
t = load(in4);

% set needed variables
z = dparam(:,1) - offset; 

% t = decaycurve(1,:); 
% S = decaycurve(2:end, :); 
%t = decaycurve(:,1); 

S = decaycurve(:,2:end); 

%Dk = olddat(:,4)*1.16e-5; 
%z_dk = olddat(:,1); 

Dk = DPPdat(:,2)*1.16e-5; % converts K from m/day to m/s
z_dk = DPPdat(:,1); 

%% Match depths of permeability measurements from my old data
% with the depths in the new processed data -- need to find nearest
% neighbors, exlude NMR points that do not have corresponding k
% measurement. 

for i = 1:length(Dk)
   [errorz(i), ind(i)] = min(abs(z_dk(i) - z)); 
end

%% Set rest of variables
phi = dparam(ind,2); % totalf
ClayH20 = dparam(ind,3); % clayf
CapH20 = dparam(ind,4); % capf
FreeH20 = dparam(ind,5); % freef
T2ML = dparam(ind,6); %mlT2
k_SDR_VC = dparam(ind,7); %Ksdr 
k_TC = dparam(ind,8); %Ktc 
k_SOE = dparam(ind,9); % Ksoe
Tsdr = dparam(ind,10); % Tsdr
Ttc = dparam(ind,11); % Ttc
Tsoe = dparam(ind,12); % Tsoe
soe = dparam(ind,13); % soe
noise = dparam(ind,14); %noise
%soe_3s = dparam(ind,9); 
%soe_twm = dparam(ind,10); 
%soe_twm_3s = dparam(ind,11); 

z = z(ind); 

%% check SOE and T2ML values from old data and new processed data. 


%%%%%%%%%%% Plot actual decay curves at each depth
    figure(1);
    hold on
    plot(t,S,'.')
    xlabel('time')
    ylabel('S(t)')
    %legend(num2str(z))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5    
%%
% create text file output for use with NMR_script_main
%datamat = [z(:), T2ML(:), phi(:), Dk(:), ...
%    soe, soe_3s(:), soe_twm(:), soe_twm_3s(:)]; 

datamat = [z(:), T2ML(:), phi(:), Dk(:), ...
    soe]; 

save([cd '/NMR_dat_processed/' name '.txt'], 'datamat', '-ascii')
end