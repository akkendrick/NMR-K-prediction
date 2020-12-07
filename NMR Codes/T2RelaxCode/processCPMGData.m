% Invert CPMG data from a data.1d file

%clear all

figNum = 5;
decay = 1;
norm = 0;

data=load('data_ZEOX.csv');
% Convert to seconds
t=1e-3*data(:,1);

% R = data(:,2)./data(1,2);
% I = data(:,3)./data(1,2);

R = data(:,2);
I = data(:,3);

P=R+1i*I;
% [a,b,c]=rotate_data(P);
% 
%  m1=mean(a(end-500:end));
%  m2=mean(b(end-500:end));
% 
% R=a-m1;
% I=b-m2;

% Create acceptable values of T2 vector (in seconds)
T2=logspace(-3,1,200);
% Set smoothing parameter
alpha=100;

RprocessCPMG = R;

if length(R)<=80000
[m,DSYN1] = whittallmodt2fit(R',t',I',T2,alpha);
else
[m,DSYN1] = whittallmodt2fit(R(1:2:end)',t(1:2:end)',I(1:2:end)',T2,alpha);    
end

if decay == 1
    figure(figNum+20)
    hold all
    plot(t,R)
    plot(t,I)
end
    

S=max(R)
maxAmp = max(m);
T2ml = 10^(sum(m(1:200).*log10(T2(1:200)'))/sum(m(1:200)))

if norm == 0 
    %m = m / maxAmp;
    m = m;
else
    m = m / norm;
end

figure(2)
hold all
semilogx(T2,m(1:200),'LineWidth', 2)
title('Calculated T_2 distribution ')
xlabel('T_2 , ms ')
xlim([min(T2),max(T2)])
set(gca,'xscale','log')
