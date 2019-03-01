% Build Wisconsin data files

names = {'G6_W2_tr5_20x_16p75_up_F_wRIN_wRFI_reg50_Va1','G5_W1_tr5_20x_16p5_up_F1n2_wRIN_wRFI_Reg50_Va1',...
    'Pl_W1_Tr5_20x_MPp75aLS_F1n2_wRIN_wRFI_Reg50_Va1','W2_Tr5_20x_MPp75aLS_Reg50_wRIN_wRFI_Va1'};
sites = {'Site1-WellG6','Site1-WellG5','Site2-WellPN1','Site2-WellPN2'};

% Get conversion factor to take NMR data from reference point to top of
% ground surface for direct comparision to DPP measurements
depthOffsets = [0.75,0.95,0.75,0.75];

for kk = 1:length(names)
    build_nmr_txt_file(sites{kk},names{kk},depthOffsets(kk))
end

