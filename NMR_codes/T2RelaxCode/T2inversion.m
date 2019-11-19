%clear all
data=load('data_distilled_CPMG.csv');
% Convert to seconds
t=1e-3*data(:,1);

%%
% Rotate real and imaginary
R=data(:,2)./data(1,2);
I=data(:,3)./data(1,2);
 P=R+1i*I; 
[a,b,c]=rotate_data(P);
 
 m1=mean(a(end-500:end));
 m2=mean(b(end-500:end));
% % 
% 
% % Subtract mean, only do this if T2 has fully decayed
R=a-m1;
I=b-m2;

%%
% background=load('background2.csv');
% %background files are all for teflon sample holders with parafilm. 1 is the
% %old background file, 2 is for a large sampleholder, from June 2015, and 3
% %is for a small sample holder, June 2015. 2 and 3 are almost identical.
% Y=1;
% 
% if background(1)==data(1)&&Y==1
%    
%     if length(data)>=length(background)
%       BG=zeros(1,length(data));
%       BG(1:length(background))=background(:,2);
%        R=R-BG';
%     else if length(data)<=length(background)
%        BG=background(1:length(data),2);
%        R=R-BG';
%         end
%     end
%     
% end


%%

figure
plot(t,R)
hold all
plot(t,I)

%%
% Create acceptable values of T2 vector (in seconds)
T2=logspace(-3,1,200);
% Set smoothing parameter
alpha=100;

if length(R)<=80000
[m,DSYN1] = whittallmodt2fit(R',t',I',T2,alpha);
else
[m,DSYN1] = whittallmodt2fit(R(1:2:end)',t(1:2:end)',I(1:2:end)',T2,alpha);    
end
%%
S=max(R)
T2ml= 10^(sum(m(1:200).*log10(T2(1:200)'))/sum(m(1:200)))
%%
save('T2spectrum.mat','T2','m','DSYN1','alpha','S','T2ml');
%%
figure(2)
semilogx(T2,m(1:200),'LineWidth', 2)
title('Calculated T_2 distribution ')
xlabel('T_2 , s ')
hold on
xlim([1e-3,10])

%saveas(gcf,'T2_dist.fig')
%saveas(gcf,'T2_dist.eps')
%%
