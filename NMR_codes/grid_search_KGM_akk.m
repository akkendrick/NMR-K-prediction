function [bestTau, bestRho, r] = grid_search_KGM(T2, phi, kk)
% this function does a grid search over specified parameter values in the
% SDR equation

if nargin < 4
    n = 2; 
end

% set bounds on rho and tau 
max_rho = 0;
min_rho = -6;

max_tau = 10; 
min_tau = 1; 

nm = 1000;
rhoSpace = logspace(min_rho, max_rho, nm); 
tauSpace = linspace(min_tau, max_tau, nm); 

% SDR equation
%SDR_f = @(b, m, n, lT2, lphi) log10(b) + m*lphi + n*lT2; 

temp = 20;
T_B = 3.3 + 0.044*(temp - 35);
eta = 1002*10^-6;
D = (1.0413+0.039828*temp+0.00040318*temp^2)*10^-9;
density = 1000*(1-(temp+288.94)/(508929*(temp+68.12))*(temp-3.98)^2);
g = 9.81; %m/s^2

KGM_f = @(rho, tau) density.*g/(8*tau^2*eta).*phi.*(-D/rho+sqrt((D/rho)^2+...
    4*D*(T2.^(-1)-T_B^(-1)).^(-1)) ).^2;

% run grid search
[r] = deal(zeros(length(tauSpace), length(rhoSpace)));  

for rholoop = 1:length(rhoSpace)
    rho = rhoSpace(rholoop);    
    for tauloop = 1:length(tauSpace)    
        tau = tauSpace(tauloop);  
        Kp = KGM_f(rho, tau);
        r(rholoop, tauloop) = norm(Kp - kk); 
    end
end

% compute best-fitting values
s = size(r); 
bestp = min(r(:)); 
indp = find(r(:) == bestp); 
[ii, jj] = ind2sub(s, indp); 

% Plot grid results
figure; 
hold on
imagesc(tauSpace, log10(rhoSpace), r./length(kk))
plot(tauSpace(jj),log10(rhoSpace(ii)),'*w','MarkerSize',10)
colorbar
colormap(jet)
ylim([log10(min(rhoSpace)), log10(max(rhoSpace))])
xlim([min(tauSpace), max(tauSpace)])
ylabel('log_{10} rho')
xlabel('tau')
%set(gca, 'YDir','reverse')

str = ['(a)  Misfit - best value pair: tau = ', num2str(tauSpace(jj)) ', log(rho) = ',num2str(log10(rhoSpace(ii)))];
title(str)


set(gcf,'OuterPosition',[100 100 1000 800], 'PaperPositionMode','auto','color',[1 1 1])

% output best b values
bestTau = tauSpace(jj)
bestRho = rhoSpace(ii)
% ms = [mspace(jj); [0,1, 2, 4]']; 
% bs = [log10(bspace(jj)); bestb(:)];
% ms = [mspace(ii); [1, 2, 4, 0]'];

end