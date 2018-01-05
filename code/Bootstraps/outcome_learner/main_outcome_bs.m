function []=main_outcome_bs()

clear; close all; clc
pati=[fileparts(fileparts(fileparts(pwd))) '/data/'];
addpath([fileparts(fileparts(fileparts(pwd))) '/code/helpfunctions/']);

%get files
dum1=strcat(pati,'Group1*.mat');    %Knockout group (Cav1.2NesCre)
dum2=strcat(pati,'Group2*.mat');    %Control (Cav1.2flfl)
files1=dir(dum1); files2=dir(dum2);
files=[files1; files2];
[nFiles, ~]=size(files);

%loop over animals
for i=1:nFiles
    
    %load data
    file=[pati files(i).name];
    disp(file)
    vpn=file(end-5:end-4);
    vpn(vpn=='_')=[];
    vpn=str2num(vpn);
 
    m=load(file);
    mtx=m.mtx(:,1:6);
    t2=m.tp(1)+1; t3=m.tp(2)+1;
    %if vpn==17; t2=104; t3=104; end
    days=1:mtx(end,1);

    gg(i)=m.g;  %group
    vps(i)=vpn; %ID
    
    %get accuracy and strategies of animal
    for iday=1:days(end)
        
        ind=mtx(:,1)==iday;
        mtxd= mtx(ind,:); T=size(mtxd,1);        
        corr=mtxd(:,4);
        response=mtxd(:,5);
            
        s_pmtx=zeros(7,T-1);
        rp_pmtx=zeros(6,T-1);
        nT(iday)=T;
        acc(iday)=sum(corr)/T;
        ind=mtxd(:,3)==1|mtxd(:,3)==2;
        accstay(iday)=sum(corr(ind))/sum(ind);
        ind=mtxd(:,3)==3|mtxd(:,3)==5;
        accshift(iday)=sum(corr(ind))/sum(ind);
        
        for t=2:T
            resp=response(t); lastresp=response(t-1);
            tmp=resp==lastresp;
            rew=corr(t); lastrew=corr(t-1);
            cond=mtxd(t,3); trial=mtxd(t,2); lasttrial=mtxd(t,2);
            
            %stay and shift strategies (for reward probabilities)
            [s_pmtx,rp_pmtx]=updateStrats(tmp, rew, lastrew,s_pmtx,rp_pmtx,t);
        end
        [~,sps]=getRP(rp_pmtx,s_pmtx);

        swinstay(iday)=sps(3);
        swinshift(iday)=sps(4);
        slosestay(iday)=sps(5);
        sloseshift(iday)=sps(6);
        
        clear s_* rp_*
    end
        
    %shown strategies of animal for simulation
    stratvals(:,1)=swinstay';
    stratvals(:,2)=swinshift';
    stratvals(:,3)=slosestay';
    stratvals(:,4)=sloseshift';
    tp=m.tp;
        
    %compare accuracy via bootstrap distribution
    nsim=1000;
    for isim=1:nsim
       [accuracy(isim,:), accuracyshift(isim,:), accuracystay(isim,:)]=sim_outcome(stratvals,tp,nT);
    end
    
     %plots actual accuracy vs. accuracy of bootstrap learner
    alpha=.05;
    down=ceil(alpha*nsim);
    up=nsim-ceil(alpha*nsim);
    labs={'accuracyshift','accuracystay'};
    labs2={'accshift',' accstay'};
    h1=figure('color','white'); hold on; box on; fs=10;
    
    for j=1:length(labs)
        eval(['dum=' labs{j} ';'])
        eval(['dumdata=' labs2{j} ';'])
        
        for iday=1:size(acc,2)
            tmp=sort(dum(:,iday),'ascend');
            ci(iday,:)=[tmp(down) tmp(up)];
        end
        maccuracy=mean(dum);
        tmpx(:,2)=(maccuracy'-ci(:,1));
        tmpx(:,1)=abs(maccuracy'-ci(:,2));
        x=1:size(accuracy,2);
        subplot(1,2,j); hold on; box on
        cols=[.98 .98 .98; .95 .95 .95;  .92 .92 .92];

        xxx=[1 tp(1) tp(2) tp(3)];
        for k=1:3
            x1=xxx(k); x2=xxx(k+1);
            h=area([x1 x2],[.995 .995],'FaceColor',cols(k,:),'EdgeColor','none');
            hpatch = get(h,'children'); set(hpatch,'FaceAlpha',0.8)
        end
        shadedErrorBar(x,maccuracy,tmpx,'b');
        plot(x,maccuracy,'b','LineWidth',2);
        plot(x,dumdata,'-r','LineWidth',1.5);
        xlabel('day#','FontSize',fs); ylabel('% correct','FontSize',fs);
        glabs={'Ca_v1.2{NesCre}','Ca_v1.2{fl/fl}'}; 
        sample=glabs{m.g};
        tit=sprintf('%s ID:%d %s',sample,vpn,labs2{j});
        title(tit,'FontSize',fs);
        set(gca,'FontSize',fs);
        xlim([0 tp(3)])
    end
    keyboard
    
    
    close
    clearvars -except nFiles files i pati crepos grouplab pat* Res* pmtx* gg vps
    
end


%--------------------------------------------------------------------------
%collect evidence for strategies, and get %of applied strategies
function [s_pmtx,rp_pmtx]=updateStrats(tmp, rew, lastrew,s_pmtx,rp_pmtx,t)
 
strats={'s_stay','s_switch','s_winstay','s_winshift','s_losestay','s_loseshift','s_corr'};
rps={'rp_stay','rp_switch','rp_winstay','rp_winshift','rp_losestay','rp_loseshift'};
        
switch tmp
    case true %same response (stay)
        s_pmtx(ismember(strats,'s_stay'),t-1)=1;
        
        if rew==1
            s_stay(t-1)=1;
            rp_pmtx(ismember(rps,'rp_stay'),t-1)=1;
        else s_stay(t-1)=-1;
            rp_pmtx(ismember(rps,'rp_stay'),t-1)=-1;
        end
        
        if lastrew==1
            s_winstay(t-1)=1;
            s_pmtx(ismember(strats,'s_winstay'),t-1)=1;
            if rew==1
                rp_pmtx(ismember(rps,'rp_winstay'),t-1)=1;
            else
                rp_pmtx(ismember(rps,'rp_winstay'),t-1)=-1;
            end
            
        elseif lastrew==0
            s_losestay(t-1)=1;
            s_pmtx(ismember(strats,'s_losestay'),t-1)=1;
            
            if rew==1
                rp_pmtx(ismember(rps,'rp_losestay'),t-1)=1;
            else
                rp_pmtx(ismember(rps,'rp_losestay'),t-1)=-1;
            end
        end
        
    case false %diff response (shift)
        s_sshift(t-1)=1;
        s_pmtx(ismember(strats,'s_shift'),t-1)=1;
        
        if rew==1,
            s_shift(t-1)=1;
            rp_pmtx(ismember(rps,'rp_shift'),t-1)=1;
        else
            s_shift(t-1)=-1;
            rp_pmtx(ismember(rps,'rp_shift'),t-1)=-1;
        end
        
        if lastrew==1 %
            s_winshift(t-1)=1;
            s_pmtx(ismember(strats,'s_winshift'),t-1)=1;
            if rew==1
                rp_pmtx(ismember(rps,'rp_winshift'),t-1)=1;
            else
                rp_pmtx(ismember(rps,'rp_winshift'),t-1)=-1;
            end
            
        elseif lastrew==0
            s_loseshift(t-1)=1; %lose-shift
            s_pmtx(ismember(strats,'s_loseshift'),t-1)=1;
            if rew==1
                rp_pmtx(ismember(rps,'rp_loseshift'),t-1)=1;
            else
                rp_pmtx(ismember(rps,'rp_loseshift'),t-1)=-1;
            end
            
        end
end


%--------------------------------------------------------------------------
%get reward probabilities/shown strategies per day
function  [rps,sps,spstot]=getRP(rp_pmtx,s_pmtx)

%burnin trials
try
    rp_pmtx(:,11:end);
    s_pmtx(:,11:end);
end

%reward probabilities
for j=1:size(rp_pmtx,1)
   tmp=rp_pmtx(j,:);
   n=sum(tmp~=0);       
   npos=sum(tmp==1);    
   rps(j)=npos/n;
end

%stay shift
ntot=size(s_pmtx,2);
nstay=sum(s_pmtx(1,:)==1);
nshift=sum(s_pmtx(2,:)==1);
n=nstay+nshift;
sps(1)=nstay/n;
sps(2)=nshift/n;

spstot(1)=nstay/ntot;
spstot(2)=nshift/ntot;

%winstay winshift
nwinstay=sum(s_pmtx(3,:)==1);
nwinshift=sum(s_pmtx(4,:)==1);
n=nwinstay+nwinshift;
sps(3)=nwinstay/n; 
sps(4)=nwinshift/n;

spstot(3)=nwinstay/ntot;
spstot(4)=nwinshift/ntot;

%losestay loseshift
nlosestay=sum(s_pmtx(5,:)==1);
nloseshift=sum(s_pmtx(6,:)==1);
n=nlosestay+nloseshift;
sps(5)=nlosestay/n;
sps(6)=nloseshift/n;

spstot(5)=nlosestay/n;
spstot(6)=nloseshift/ntot;

ind=sps==0;
rps(ind)=0;
if sum(isnan(rps))>0,keyboard; end


