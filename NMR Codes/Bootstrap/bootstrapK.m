function [bootstat] = bootstrapK(NBoot,K_DPP,K_NMR)
%BOOTSTRAPK Bootstrap K estimates to calculate the best SDR based b value.
%The function assumes K_NMR is a SDR derived model like K parallel or K
%series containing a K SDR with a b value somewhere

disp('Test')

data = [K_DPP, K_NMR']
%bootstat = bootstrp(NBoot,min(K_DPP-K_NMR'),K_NMR');
bootstrp(NBoot,@(x) bootstrap_NMR_K(x), data);
disp('Test')


end

