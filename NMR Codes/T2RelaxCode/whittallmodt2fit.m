function [m,dsyn] = whittallmodt2fit(d,t,noise,T2,eps)

%whittallmodt2fit.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----Inversion of NMR data for T2 distribution (smoothness)----%
%Written by                                                    %
%Elliot Grunewald                                              %
%September 2006                                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                                             
% Inverts time series NMR data for a distribution of multi exponential
% decays.  The problem is linearized by predetermining a series of T2
% values to which the distribution can be fit.  Data are resampled to a log
% spacing using resampled vectors rd and rt.  In the inversion, data and
% corresponding rows of the L operator are normalized by the standard
% deviation of the corresponding datum.  The inversion is also regularized
% by a smoothness constraint. DATA MAY BE FULLY DECAYED OR NOT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INPUT ARGUMENTS
%[m,dsyn] = whittallmodt2fit(d,t,noise,T2,eps)
%d     data vector corresponding to t  !!!!!!!NO TIME ZERO ALLOWED!!!!!!!!
%t     time vector corresponding to d 
%noise noise vector corresoponding to t and d
%T2    vector of available T2 values to which the distribution is fit
%eps   regularization factor 
%
%%OUTPUT ARGUMENTS
%m     inverted model 1,201 with m(201)=baseline offset
%dsyn  modeled data with dysn(:,1)=resampled t and dsyn(:,2)=M(t)


%%Get data statistics
stdev=std(noise); %%set stdev per point to stdev of channel

%NOTE THAT THE BELOW LINE IS NOT USUALLY COMMENTED, DENYS COMMENTED THIS
%AUG 22, 2011

disp(['Signal-to-Noise = ' num2str(d(1)/stdev)]); %%display signal to
%noise ratio        

%n=20;  %%down sampling
%n=5;
n=1;

%%%%%%%%%%%%%%%%
%NOTE THIS GOT CHANGED, PROBABLY ERROR HERE IF ERROR
if length(d)<1500
    n=1;
end
%%%%%%%%%%%%%


disp(['down sampling by factor of ' num2str(n)]); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%NOTE THIS CAN BE TURNED OFF OR ON IT STILL WORKS IF ITS GIVING YOU AN ERROR
%ITS PROBABLY HERE

% %clip decayed data
% %by finding when moving average reaches some fraction of initial signal
% winsize=20;  %%define size of moving average window
% filt=1/winsize*ones(1,winsize);%%create moving average filter
% %enddecay=min(find(conv(filt,d)<d(1)/1000));%%find when mov avg < 1/1000 of d(1) 
% enddecay=min(find(conv(filt,d)<d(1)/500)); %%use for MIXCPMG
% d=d(1:enddecay);%%trim d
% t=t(1:enddecay);%%trim t
% noise=noise(1:enddecay);%%trim noise
% % UNTIL HERE, LEAVE EVERYTHING ELSE ALONE

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Reduce sampling
disp('Resampling data')
tlog_ideal=logspace(log10(t(length(t))),log10(t(1)), length(t)/n);  %the ideal log spacing of data
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
        BB=stdev/sqrt(cindex-lastindex);
        tempstdev(i)=stdev/sqrt(cindex-lastindex); %stdev is std divided by number of data in the new binned point
    end
    lastindex=cindex;
end

nonnan=(find(abs(tempt)>0));  %%clear out all points in resampled data which are NaN
rt=tempt(nonnan)';
rd=tempd(nonnan)';
resampleddata=rd;
rstdev=tempstdev(nonnan)';
rd=rd./rstdev;%%weight the data by their standard deviation

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
%Lwoweightreg=[Lwoweightreg, ones(length(rt),1)];
%L=[L,ones(length(rt),1)]; %%add column of ones for baseline offset

%%Add regularization rows
k=length(T2)+1; %%number of columns in L
reg = diag(-2*ones(k,1),0)+diag(ones(k-1,1),1)+diag(ones(k-1,1),-1);
reg = [reg(2:k-2,1:k-1),zeros(k-3,1)];
Lreg = [L;(eps*reg)]; %%append the regularization matrix to L

%%Calcuate inversion using non-negative lsq
disp('inverting data')
interm_data=Lreg;
[m]=lsqnonneg(Lreg,[rd;zeros(k-3,1)]); %append zeros to d %%SLOWER
%options=optimset('TolFun',1e-25);
%[m]=lsqlin(Lreg,[rd;zeros(k-3,1)],[],[],[],[],zeros(201,1),[],[],options);  %%FASTER
dsyn(:,1)=rt;
dsyn(:,2)=Lwoweightreg*m;
dsyn(:,3)=resampleddata-dsyn(:,2);
%%WRONG!!disp(['Data tolerance =' num2str(stdev*sqrt(length(rd)))])
disp(['Data tolerance =' num2str(sqrt(sum(rstdev.^2)))])
disp(['Done -- Norm Residual = ' num2str(norm(dsyn(:,3)))])
%residual=resampleddata-Lwoweightreg*m;

