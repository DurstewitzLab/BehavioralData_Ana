function []=main_cue_bs()
%start cue bootstraps

clear; close all; clc
pati=[fileparts(fileparts(fileparts(pwd))) '/data/'];
addpath([fileparts(fileparts(fileparts(pwd))) '/code/helpfunctions/']);

%get files
dum1=strcat(pati,'Group1*.mat');    %Knockout group (Cav1.2NesCre)
dum2=strcat(pati,'Group2*.mat');    %Control (Cav1.2flfl)
files1=dir(dum1); files2=dir(dum2);
files=[files1; files2];
[nFiles, ~]=size(files);

for i=1:nFiles
    file=[pati files(i).name];
    disp(file)
    vpn=file(end-5:end-4);
    vpn(vpn=='_')=[];
    vpn=str2num(vpn);
 
    m=load(file);
    mtx=m.mtx(:,1:6);
    t2=m.tp(1)+1; t3=m.tp(2)+1;
    if vpn==17; t2=104; t3=104; end
    days=1:mtx(end,1);

    gg(i)=m.g;  %group
    vps(i)=vpn; %ID
    
    %get reward probabilities, strategies, and accuracies per day
    for iday=1:days(end)
        
        ind=mtx(:,1)==iday;
        mtxd= mtx(ind,:); 
      
        T=size(mtxd,1);
        corr=mtxd(:,4);
        response=mtxd(:,5);     %left right response
        cues=mtxd(:,6);         %top bottom cue
        conds=mtxd(:,3);
        
        s_pmtx=zeros(4,T);
        rp_pmtx=zeros(4,T);
        sd_pmtx=zeros(4,T);
        rpd_pmtx=zeros(4,T);
        
        %get accuracy for comparison
        nT(iday)=T;
        acc(iday)=sum(corr)/T;
        ind=mtxd(:,3)==1|mtxd(:,3)==2;
        accstay(iday)=sum(corr(ind))/sum(ind);
        ind=mtxd(:,3)==3|mtxd(:,3)==5;
        accshift(iday)=sum(corr(ind))/sum(ind);
        
        for t=1:T
            resp=response(t); cue=cues(t); rew=corr(t); cond=conds(t);
            if cond==4, cue=3; end
            
            %stay and shift strategies (for reward probabilities)
            if cond==1 || cond==3 
                [s_pmtx,rp_pmtx]=updateStrats(resp,cue,rew,s_pmtx,rp_pmtx,t);
            elseif cond==2 || cond==5 || cond==4
                [sd_pmtx,rpd_pmtx]=updateStratsD(resp,cue,rew,sd_pmtx,rpd_pmtx,t);
            end
        end
        [rps,sps]=getRP(rp_pmtx,s_pmtx);
        [rpsd, spsd]=getRP(rpd_pmtx, sd_pmtx);
        
        rpss(iday,:)=rps;
        spss(iday,:)=sps;
        rpssd(iday,:)=rpsd;
        spssd(iday,:)=spsd;     
        
        clear s_* rp_* sd_* rpd_*
    end

    tp=m.tp; xx=[1 109 124];
    
    stratvals1=spss; tp=m.tp;
    %if vpn==17, tp=[105 103 103]; end
    stratvals2=spssd;

    %performs bootstraps and returns accuracy of simulated learner
    nsim=1000;
    for isim=1:nsim
        [accuracy(isim,:), accuracyshift(isim,:), accuracystay(isim,:)]=sim_cue(stratvals1,stratvals2,tp,nT);
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
    clearvars -except nFiles files i pati crepos pat* Res* pmtx* gg vps
    
end




%--------------------------------------------------------------------------
%collect evidence for strategies, and get %of applied strategies
function [s_pmtx,rp_pmtx]=updateStrats(resp,cue, rew, s_pmtx,rp_pmtx,t)

%1: cue=1, resp=1   %incorrect
%2: cue=1, resp=0   %c
%3: cue=2, resp=1   %c
%4: cue=2, resp=0  %inccorrect

switch resp
    case 0
        if cue==1
            s_pmtx(2,t)=1;
            
            if rew==1
                rp_pmtx(2,t)=1;
            else
                rp_pmtx(2,t)=-1;
            end
            
        elseif cue==2
            s_pmtx(4,t)=1;
            if rew==1
                rp_pmtx(4,t)=1;
            else
                rp_pmtx(4,t)=-1;
            end
        end
        
    case 1
        if cue==1
            s_pmtx(1,t)=1;
            
            if rew==1
                rp_pmtx(1,t)=1;
            else
                rp_pmtx(1,t)=-1;
            end
            
        elseif cue==2
            s_pmtx(3,t)=1;
            if rew==1
                rp_pmtx(3,t)=1;
            else
                rp_pmtx(3,t)=-1;
            end
        end
end


%--------------------------------------------------------------------------
%collect evidence for strategies, and get % applied strategies
function [s_pmtx,rp_pmtx]=updateStratsD(resp,cue, rew, s_pmtx,rp_pmtx,t)

%1: cue=1, resp=1   %incorrect
%2: cue=1, resp=0   %c
%3: cue=2, resp=1   %c
%4: cue=2, resp=0  %inccorrect

switch resp
    case 0
        if cue==1
            s_pmtx(2,t)=1;
            
            if rew==1
                rp_pmtx(2,t)=1;
            else
                rp_pmtx(2,t)=-1;
            end
            
        elseif cue==2
            s_pmtx(4,t)=1;
            if rew==1
                rp_pmtx(4,t)=1;
            else
                rp_pmtx(4,t)=-1;
            end
        end
        
    case 1
        if cue==1
            s_pmtx(1,t)=1;
            
            if rew==1
                rp_pmtx(1,t)=1;
            else
                rp_pmtx(1,t)=-1;
            end
            
        elseif cue==2
            s_pmtx(3,t)=1;
            if rew==1
                rp_pmtx(3,t)=1;
            else
                rp_pmtx(3,t)=-1;
            end
        end
end




%--------------------------------------------------------------------------
%get reward probabilities/shown strategies per day
function  [rps,sps]=getRP(rp_pmtx,s_pmtx)

%1: cue=1, resp=1   %incorrect
%2: cue=1, resp=0   %c
%3: cue=2, resp=1   %c
%4: cue=2, resp=0  %inccorrect

%reward probabilities
for j=1:size(rp_pmtx,1)
    tmp=rp_pmtx(j,:);
    n=sum(tmp~=0);       %von allen ausgeführten
    npos=sum(tmp==1);    %wieviele davon waren successful
    rps(j)=npos/n;
end

%no of left/right responses to left/right cues
ntot=size(s_pmtx,2);
ncue1resp1=sum(s_pmtx(1,:));
ncue1resp0=sum(s_pmtx(2,:));
ncue2resp1=sum(s_pmtx(3,:));
ncue2resp0=sum(s_pmtx(4,:));
ncue1=sum(s_pmtx(1,:)+s_pmtx(2,:));
ncue2=sum(s_pmtx(3,:)+s_pmtx(4,:));

sps(1)=ncue1resp1/ncue1;
sps(2)=ncue1resp0/ncue1;
sps(3)=ncue2resp1/ncue2;
sps(4)=ncue2resp0/ncue2;

ind=sps==0;
rps(ind)=0;
