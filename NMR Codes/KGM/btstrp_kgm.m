function [out] = btstrp_kgm(X,m)
% this function does a grid search over specified parameter values in the
% SDR equation

lT2 = X(:,1);
lphi = X(:,2);
kk = X(:,3);

% put in linear space
T2 = 10.^lT2; 
phi = lphi; 
K = kk; 

% temp-dependent parameters
Tb= @(T) 3.3 + .044.*(T - 35);                  % Bulk relaxation time constant (s)
eta = @(T) 2.414e-5*(10.^((247.8)./(T + 133.15)));  % viscosity (Pa s)
aa = 1.0413; bb = 0.039828; cc = 0.00040318;            % diffusivity constants
D = @(T) aa + bb*T+ cc*T.^2;   % diffusivity (m^2/s)
rho_h2o = @(T) 1000*(1 - (((T + 288.94)/(508929 ...
    *(T+68.12)))*((T - 3.98).^2)));             % density of water (kg/m^3)

% set bounds on A = g/(8*tau^2*rho^2), B = rho^2
max_tau = 0; 
min_tau = -2; 

max_rho = 2; 
min_rho = -6.6; 

nm = 201; 
tauspace = linspace(min_tau, max_tau, nm); 
rhospace = logspace(min_rho, max_rho, nm); 

parameterMatrixTau = tauspace .* ones(nm,nm);
parameterMatrixRho = rhospace' .* ones(nm,nm);

parameterMatrixTau = (fliplr(parameterMatrixTau));
parameterMatrixRho = (fliplr(parameterMatrixRho));

% KGM equation
% lK_kgm = @(A,B, T, T2, lphi) (A) + lphi ...
%     + 2*log10( - (D(T)/B) + sqrt((D(T)/B)^2 + ((4*D(T)*Tb(T)*T2)./(Tb(T) - T2)))); 
%lK_kgm = @(A, B, T, T2, phi) KGMfun([A, B], [T2]); 
lK_kgm = @(A, B, C, T, T2, phi) KGMfun([A, B, C], [T2, phi]); 

T = 11.111; 
% run grid search
[r] = deal(zeros(length(tauspace), length(rhospace)));  
for mloop = 1:length(rhospace)
    for bloop = 1:length(tauspace)   
        btest = parameterMatrixRho(mloop,bloop);
        atest = parameterMatrixTau(mloop,bloop);  
        lKp = lK_kgm(atest, btest, m, T, T2, phi);
        r(mloop, bloop) = norm(lKp - kk); 
    end
end

% % compute best-fitting values
% ind1 = Bspace == 1; 
% ind2 = Bspace == 2; 
% ind4 = Bspace == 4; 
% ind0 = Bspace == 0; 
% 
% r1 = r(ind1, :) == min(r(ind1,:)); 
% r2 = r(ind2, :) == min(r(ind2,:)); 
% r4 = r(ind4,:) == min(r(ind4,:)); 
% r0 = r(ind0,:) == min(r(ind0,:)); 
% 
% bestb(1)= log10(Aspace(r1)); 
% bestb(2)= log10(Aspace(r2)); 
% bestb(3)= log10(Aspace(r4)); 
% bestb(4)= log10(Aspace(r0)); 

s = size(r); 
bestp = min(r(:)); 

indp = find(r(:) == bestp); 
[ii, jj] = ind2sub(s, indp);
% bestA = tauspace(ii); 
% bestB = rhospace(jj); 

bestA = parameterMatrixTau(ii,jj);
bestB = parameterMatrixRho(ii,jj);

lkpred = lK_kgm(bestA, bestB, m, T, T2, phi); 

taus = 1./sqrt(10.^(tauspace)); 
tausParamMatrix = 1./sqrt(10.^parameterMatrixTau);
rhosParamMatrix = log10(parameterMatrixRho);

errorMatrix = log10(r/length(phi));
% display('ii')
% ii
% display('jj')
% jj
% 
% % Make sure all variables are displayed in LINEAR SPACE
% display('tau (ii,jj)')
besttau = tausParamMatrix(ii,jj);
% display('rho (ii,jj)')
bestrho = bestB;
% display('Parameter Error (ii,jj)')
% errorMatrix(ii,jj)

% KGM_lk = lkpred;
% 
% % tausParamMatrix = fliplr(tausParamMatrix);
% % rhosParamMatrix = flipud(rhosParamMatrix);
% % errorMatrix = flipud(fliplr(errorMatrix));
% bestp = min(errorMatrix(:)); 
% 
% indp = find(errorMatrix(:) == bestp); 
% [ii, jj] = ind2sub(s, indp);


KGM_lk = lkpred;

% % Plot grid results
% figure; 
% % subplot(4, 4, [1:3,5:7,9:11, 13:15])
% hold on
% %imagesc(fliplr(taus), fliplr(log10(rhospace)), errorMatrix)
% %imagesc([1 10],[0 -6.6],flipud(fliplr(errorMatrix)))
% %imagesc([1 10],[0 -6.6],errorMatrix)
% 
% %uimagesc(fliplr(taus),fliplr(log10(rhospace)),errorMatrix)
% uimagesc(fliplr(taus),(log10(rhospace)),errorMatrix)
% box on
% grid on
% 
% % tausParamMatrix = fliplr(tausParamMatrix);
% % rhosParamMatrix = flipud(rhosParamMatrix);
% % errorMatrix = flipud(fliplr(errorMatrix));
% 
% bestp = min(errorMatrix(:)); 
% 
% indp = find(errorMatrix(:) == bestp); 
% [ii, jj] = ind2sub(s, indp);
% 
% plot(tausParamMatrix(ii,jj), rhosParamMatrix(ii,jj), '*w')
% caxis([-1.2,0])
% colorbar('Ticks', [-1.2:.2:0])
% colormap(jet)
% xlim([min((taus)), max((taus))])
% ylim([min(log10(rhospace)), max(log10(rhospace))])
% ylabel('log_{10} (\rho)')
% xlabel('\tau')
% set(gca,'FontSize',14)


out = [besttau; bestrho];
end