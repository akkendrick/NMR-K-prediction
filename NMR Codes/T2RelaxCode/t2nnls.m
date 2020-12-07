function [m,dsyn,r,ChiSq] = t2nnls(d,t,noise,T2,eps,NSamp)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ----Inversion of NMR data for T2 distribution (smoothness)---- % 
% Written by                                                     %
% Elliot Grunewald                                               %
% September 2006                                                 %   
% REVISION 12/3/2009                                             %   
% Keating and Falzone Revisions 08/31/2012                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                                             
% Inverts time series NMR data for a distribution of multiexponential
% decays.  The problem is linearized by predetermining a series of T2
% values to which the distribution can be fit. Data are resampled to a log
% spacing using resampled vectors rd and rt.  In the inversion, data and
% corresponding rows of the L operator are normalized by the standard
% deviation of the corresponding datum.  The inversion is also regularized
% by a smoothness constraint. DATA MAY BE FULLY DECAYED OR NOT
% 08/31/2012 revisions
%  -- updated regularization matrix                              %
%  -- added moving average data filter    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INPUT ARGUMENTS
%[m,dsyn,r ] = t2nnls(d,t,noise,T2,eps)
%d     data vector corresponding to t  !!!!!!!NO TIME ZERO ALLOWED!!!!!!!!
%t     time vector corresponding to d 
%noise noise vector sqcorresoponding to t and d
%T2    vector of available T2 values to which the distribution is fit
%eps   regularization factor 
%NSamp down sampling #
%
%%OUTPUT ARGUMENTS
%m     inverted model 1, length(T2)+1 with m(length(T2)+1)=baseline offset
%dsyn  modeled data with dysn(:,1)=resampled t and dsyn(:,2)=M(t) ,
%      dsyn(:,4) = ressampled d
%r      residual norm and data tol r(1) = residual norm; r(2) = data tol

%%Get data statistics
stdev=std(noise); %%set stdev per point to stdev of channel

%disp(['Signal-to-Noise = ' num2str(d(1)/stdev)]); %%display signal to noise ratio
n=NSamp;  %%down sampling - default 20 for lab data

%%clip decayed data
%%by finding when moving average reaches some fraction of initial signal
winsize=100;  %%define size of moving average window
filt=1/winsize*ones(1,winsize);%%create movinng average filter
enddecay=min(find(conv(filt,d)<d(1)/1e5));%%find when mov avg < 1e-5 of d(1)
% if(length(enddecay)<1 || enddecay > length(d))
if(length(enddecay)<=1 || enddecay > length(d)) % modified for noisy data where d(1) is negative
    enddecay=length(d);
end

% %for test plot
% d2 =d; t2 = t; noise2 = noise;
% %for test plot

d=d(1:enddecay);%%trim d
t=t(1:enddecay);%%trim t
noise=noise(1:enddecay);%%trim noise

% %test figure to check how much of the dat aand noise are clipped
% figure('units','normalized','position',[0 0.65 .2 .25])
% subplot(211); plot(t2,d2); hold on; plot(t,d,'r'); grid; legend('Ini','Trimed')
% subplot(212); plot(t2,noise2); hold on; plot(t,noise,'r'); grid;  legend('Ini','Trimed')
% %test figure to check how much of the dat aand noise are clipped

% apply moving average filter with window size of 5
winsize = 5;
dtemp = d(end:-1:1);
dtemp = filter(ones(1,winsize)/(winsize),1,dtemp);
d=dtemp(end:-1:winsize);
t = t(1:end-winsize+1);
noise = noise(1:end-winsize+1);

% %%Reduce sampling (log sampling)
%disp('Resampling data')

tlog_ideal=logspace(log10(t(length(t))),log10(t(1)), length(t)/n);  %the ideal log spacing of data
size(tlog_ideal)
tlog_ideal=tlog_ideal(length(tlog_ideal):-1:1); %reverse order so goes from small to large
tempt=zeros(1,length(tlog_ideal)); 
tempd=zeros(1,length(tlog_ideal)); 
tempstdev=zeros(1,length(tlog_ideal)); 
tempd=d(1);
tempstdev=stdev;
tempt(1)=t(1); 
tempd(1)=d(1);
tempstdev(1)=stdev;

lastindex=1;

for i=2:length(tlog_ideal) %for each ideal time point
    cindex=find(abs(t-tlog_ideal(i))==min(abs(t-tlog_ideal(i))));  %find closest datum to ideal point
    if(cindex==lastindex) %if the closest point is the same as before assign NaN
        tempt(i)=NaN;
        tempd(i)=NaN;
        tempstdev(i)=NaN;
    else %if closest point is unique
        tempt(i)=mean(t(lastindex+1:cindex)); %resampled binned = mean for pts between last datum and closest datum
        tempd(i)=mean(d(lastindex+1:cindex));
        tempstdev(i)=stdev/sqrt(cindex-lastindex); %stdev is std divided by number of data in the new binned point
    end
    lastindex=cindex;
end

nonnan=(find(abs(tempt)>0));  %%clear out all points in resampled data which are NaN
rt=tempt(nonnan)';
rd=tempd(nonnan)';
resampleddata=rd;
rstdev=tempstdev(nonnan)';
rd=rd./rstdev;%%weight the data by their standard deviation, AB: this normalizes all data by same STD, could be optimized!

%%Setup L matrix
L=zeros(length(rt), length(T2)+1);
Lwoweightreg=zeros(length(rt), length(T2)+1);
for i=1:length(rt)
    for j=1:length(T2)
        Lwoweightreg(i,j)=exp(-rt(i)/T2(j)); %% save unweighted matrix for later use
        L(i,j)=exp(-rt(i)/T2(j))/rstdev(i); %%normalize row by data stdev
    end
        Lwoweightreg(i,length(T2)+1)=1;
        L(i,length(T2)+1)=1/rstdev(i);
end

% Created regularization matrix
k=length(T2)+1; %%number of columns in L
reg = diag(-2*ones(k,1),0)+diag(ones(k-1,1),1)+diag(ones(k-1,1),-1);
reg = [reg(2:k-2,1:k-1),zeros(k-3,1)];
%AB, CHECK THE LINE BELOW
reg = [1 zeros(1,k-1); -2 1 zeros(1,k-2); reg; zeros(1,k-3) 1 -2 1];  % changes reg matrix to [1 zeros, -2 1 zeros, 1 -2 1 zeros ... zeros 1 -2 1] -sam, 8/18/2012

Lreg = [L;(eps*reg)]; %%append the regularization matrix to L

%%Calcuate inversion using non-negative lsq

% Alex K - setting tolerance higher for DT2 processing?
options = optimset('TolX',1e-2);
[m]=lsqnonneg(Lreg,[rd;zeros(k,1)],options); %append zeros to d %%SLOWER

%options=optimset('TolFun',1e-5);                                           %%FASTER
%mi=zeros(201,1);%initial model
%mi(100)=rd(1);%all energy at middle of distribution
%[m]=lsqlin(Lreg,[rd;zeros(k-3,1)],[],[],[],[],zeros(201,1),[],[],options);  %%FASTER

dsyn(:,1)=rt;
dsyn(:,2)=Lwoweightreg*m;
dsyn(:,3)=resampleddata-dsyn(:,2);
dsyn(:,4) = resampleddata; % Raw (resampled) data for datafit plot
dsyn(:,5) = rms(resampleddata)-rms(dsyn(:,2)); % Relative RMS error
r(2) = sqrt(sum(rstdev.^2));
r(1) = norm(dsyn(:,3));
% r(3) = abs((rms(dsyn(:,2))-rms(dsyn(:,4)))/rms(dsyn(:,2)))*100;
% ChiSq added - Ahmad Behroozmand Feb 2016
ChiSq = sum((resampleddata-dsyn(:,2)).^2 ./rstdev.^2)/(length(resampleddata)-1);


% %test
% figure
% subplot(211)
% semilogx(dsyn(:,1),dsyn(:,2))
% hold on
% semilogx(dsyn(:,1),resampleddata,'k')
% grid
% subplot(212)
% plot(dsyn(:,1),dsyn(:,2))
% hold on
% plot(dsyn(:,1),resampleddata,'k')
% grid
% grid

