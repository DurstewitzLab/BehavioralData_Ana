function make_figs_2A_2B_3A_3B()



clear; close all; clc
cd ..
pati=[pwd '/data/']
addpath([pwd '/code/helpfunctions/']);


dum1=strcat(pati,'Group1*.mat');    %Knockout group (Cav1.2NesCre)
dum2=strcat(pati,'Group2*.mat');    %Control (Cav1.2flfl)
files1=dir(dum1); files2=dir(dum2);
files=[files1; files2];
[nFiles, ~]=size(files);

grouplabs={'Group1','Group2'}; %KO vs. Control

%load in all animal files and get accuracy for different event types
%--------------------------------------------------------------------------
for i=1:nFiles
    file=[pati files(i).name];
    disp(file)
    f=files(i).name; 
    vpn=f(end-5:end-4);
    vpn(vpn=='_')=[];
    vpn=str2num(vpn);
    
    %load file, and check for errors
    dat=load([pati files(i).name]);
    %%%dat.mtx explanation
    %row1:day
    %row2:trials (e.g. 100 per day)
    %row3:condition: 1: 2: 3: 4: 5
        %1: stay-HC 
        %2: stay-LC
        %3: shift-HC
        %4: ambiguous (two same light cues)
        %5: shift-LC 
    %row4: correct yes/no
    %row5:left=1, right=0
    %row6: cue for left=1, cue for right=2
    
    %----------------------------------------------------------------------
    %get information on correct trials etc
    
    mtx=dat.mtx;
    days=unique(mtx(:,1));      %no days
    res(i).rat=vpn;
    res(i).group=dat.g;
    res(i).grouplab= grouplabs{dat.g};
    res(i).phasestart=dat.tp;
    res(i).mtx=dat.mtx;
    
    %go through all days and find number of correct events vs. total events
    for iday=1:length(days)
        ind=mtx(:,1)==iday;
        mtxd=mtx(ind,:);
      
        %get accuracy (no correct vs. total number) for different events:
        labs={'StayHC','StayLC','ShiftHC','Ambig','ShiftLC','HC','LC','Shift','Stay','all'};         
        
        %accuracy for StayHC StayLC ShiftHC ShiftLC Ambiguous 
        for icond=1:5
            ind=mtxd(:,3)==icond;
            if (icond==2 || icond==4 || icond==5) && iday<2 %should not happen
                c=0; n=0;
            else
                c=sum(mtxd(ind,4));
                n=sum(ind);
            end
            eval(['res(i).' labs{icond} '_c(iday)=c;'])
            eval(['res(i).' labs{icond} '_n(iday)=n;'])
        end
        
        %accuracy only HC 
        icond=6;
        ind=mtxd(:,3)==1 | mtxd(:,3)==3;
        c=sum(mtxd(ind,4));
        n=sum(ind);
        res(i).ND_c(iday)=res(i).StayHC_c(iday)+res(i).ShiftHC_c(iday);
        res(i).ND_n(iday)=res(i).StayHC_n(iday)+res(i).ShiftHC_n(iday);
        
        %accuracy only LC
        icond=7;
        ind=mtxd(:,3)==2 | mtxd(:,3)==5;
        c=sum(mtxd(ind,4));
        n=size(mtxd(ind,:),1);
        if iday<2,
            c=0; n=0;
        else
            c=sum(mtxd(ind,4));
            n=size(mtxd(ind,:),1);
        end
        res(i).D_c(iday)=res(i).StayLC_c(iday)+res(i).ShiftLC_c(iday);
        res(i).D_n(iday)=res(i).StayLC_n(iday)+res(i).ShiftLC_n(iday);
        
        
        %all shoift trials 
        icond=8;
        ind=mtxd(:,3)==3 | mtxd(:,3)==5;
        c=sum(mtxd(ind,4));
        n=size(mtxd(ind,:),1);       
        eval(['res(i).' labs{icond} '_c(iday)=c;'])
        eval(['res(i).' labs{icond} '_n(iday)=n;'])
        
        %all stay trials 
        icond=9;
        ind=mtxd(:,3)==2 | mtxd(:,3)==1;
        c=sum(mtxd(ind,4));
        n=size(mtxd(ind,:),1);
        eval(['res(i).' labs{icond} '_c(iday)=c;'])
        eval(['res(i).' labs{icond} '_n(iday)=n;'])
        
        %all trials (all)
        icond=10;
        ind=mtxd(:,3)~=4;
        c=sum(mtxd(ind,4));
        n=size(mtxd(ind,:),1);
        eval(['res(i).' labs{icond} '_c(iday)=c;'])
        eval(['res(i).' labs{icond} '_n(iday)=n;'])
        
        clearvars -except dat nFiles files pat* res* d i* mtx vpn g  days grouplabs
    end
end

clearvars -except pat* res
%plot accuracy for KO vs. Control on selected events
createfig(res);



%--------------------------------------------------------------------------
function createfig(data)

xx=[1 109 124 143];
n=size(data,2);
labs={'StayHC','StayLC','ShiftHC','Ambig','ShiftLC','HC','LC','Shift','Stay','all'};         

%plot accuracy over the entire time course
%--------------------------------------------------------------------------
res=nan(24,2+3*10);
for opti=[8:10]   %specify trial type for figure (see labs above)
    opt=opti;
    mtx=nan(24,143);
    
    for i=1:n
        m=opti*3;
        tmp=data(i);
        g(i)=tmp.group;
        if tmp.rat==17,     %correct for rat which died after first phase
            y1=1:tmp.phasestart(1); y2=[]; y3=[]; 
        else
            y1=1:tmp.phasestart(1)-1;
            y2=tmp.phasestart(1):tmp.phasestart(2)-1;
            y3=tmp.phasestart(2):tmp.phasestart(3);
        end
        eval(['dum1=tmp.' labs{opt} '_c(y1)./tmp.' labs{opt} '_n(y1);'])
        eval(['dum2=tmp.' labs{opt} '_c(y2)./tmp.' labs{opt} '_n(y2);'])
        eval(['dum3=tmp.' labs{opt} '_c(y3)./tmp.' labs{opt} '_n(y3);'])
        
        x1=y1; x2=xx(2):xx(2)+length(y2)-1; x3=xx(3):xx(3)+length(y3)-1;
        mtx(i,x1)=dum1; mtx(i,x2)=dum2; mtx(i,x3)=dum3;
        res(i,1:2)=[tmp.rat tmp.group];
        res(i,m:m+2)=[nanmean(dum1) nanmean(dum2) nanmean(dum3)];
        
        clear dum* y*
    end
    
    %----------------------------------------------------------------------
    cols=[.98 .98 .98; .95 .95 .95;  .92 .92 .92];

    h1=figure('color','white'); hold on; fs=12; box on
    h=colormap('jet'); c1=h(5,:); c2=h(end-6,:);
    Contr=mtx(g==2,:); nContr=sum(~isnan(Contr)); mContr=nanmean(Contr); seContr=nanstd(Contr)./sqrt(nContr);
    KO=mtx(g==1,:); nKO=sum(~isnan(KO)); mKO=nanmean(KO); seKO=nanstd(KO)./sqrt(nKO);
   

    for iphase=1:3
        xarea=xx(iphase):xx(iphase+1);
        x=xx(iphase):xx(iphase+1)-1;
        h=area([xarea(1) xarea(end)],[.998 .998],'FaceColor',cols(iphase,:),'EdgeColor','none');
        hpatch = get(h,'children'); set(hpatch,'FaceAlpha',0.8)
        
        xin=~isnan(mContr(x));
        if sum(xin)>0
            l1=shadedErrorBar(x(xin),mContr(x(xin)),seContr(x(xin)),{'color',c1,'LineWidth',2})
            
            xin=~isnan(mKO(x));
            l2=shadedErrorBar(x(xin),mKO(x(xin)),seKO(x(xin)),{'color',c2,'LineWidth',2})
        end
        
    end
    set(gca,'FontSize',fs,'YTick',[0:.2:1],'YTickLabel',[0 20 40 60 80 100]); xlim([0.5 143.5])
    plot(xx(1):xx(end),ones(1,xx(end))*.5,'-.k','LineWidth',2);
    ylabel('% correct responses','FontSize',fs)
    xlabel('day #', 'FontSize',fs)
    title(['Accuracy in ' labs{opt} ], 'FontSize',fs);
    legend([l1.mainLine l2.mainLine],{'Ca_v1.2^{fl/fl}','Ca_v1.2^{NesCre}'},'FontSize',fs,'box','off')
        
    clear Contr KO nContr nKO seContr seKO mtx tmp g
end


%plot accuracy averaged over the 3 phases
%--------------------------------------------------------------------------
clear res
for opti=[8:10]   %select type of trial
    opt=opti;
    mtx=nan(24,143);

    for i=1:n
        m=opti*3;
        tmp=data(i);
        g(i)=tmp.group;
        if tmp.rat==17,
            y1=1:tmp.phasestart(1); y2=[]; y3=[]; 
        else
            y1=1:tmp.phasestart(1)-1;
            y2=tmp.phasestart(1):tmp.phasestart(2)-1;
            y3=tmp.phasestart(2):tmp.phasestart(3);
        end
        eval(['dum1=tmp.' labs{opt} '_c(y1)./tmp.' labs{opt} '_n(y1);'])
        eval(['dum2=tmp.' labs{opt} '_c(y2)./tmp.' labs{opt} '_n(y2);'])
        eval(['dum3=tmp.' labs{opt} '_c(y3)./tmp.' labs{opt} '_n(y3);'])
        
        x1=y1; x2=xx(2):xx(2)+length(y2)-1; x3=xx(3):xx(3)+length(y3)-1;
        mtx(i,x1)=dum1; mtx(i,x2)=dum2; mtx(i,x3)=dum3;
        res(i,1:5)=[tmp.rat tmp.group nanmean(dum1) nanmean(dum2) nanmean(dum3)];        
        clear dum* y*
    end
    
    %create the figure
    cols=[.98 .98 .98; .95 .95 .95;  .92 .92 .92];
    h1=figure('color','white'); hold on; fs=12; box on

    h=colormap('jet'); c1=h(5,:); c2=h(end-6,:);
    Contr=res(res(:,2)==2,3:5); nContr=sum(~isnan(Contr)); mContr=nanmean(Contr); seContr=nanstd(Contr)./sqrt(nContr);
    KO=res(res(:,2)==1,3:5); nKO=sum(~isnan(KO)); mKO=nanmean(KO); seKO=nanstd(KO)./sqrt(nKO);
        
    
    
    errorbar(1:3,mContr,seContr,'color',c1,'LineWidth',2)
    l1=plot(1:3,mContr,'color',c1,'LineWidth',2)
    errorbar(1:3,mKO,seKO,'color',c2,'LineWidth',2)
    l2=plot(1:3,mKO,'color',c2,'LineWidth',2)
    
    set(gca,'FontSize',fs,'YTick',[0:.2:1],'YTickLabel',[0 20 40 60 80 100],'XTick',1:3,'XTickLabel',{'I','II','II'}); 
    xlim([0.8 3.2]); ylim([0 1])
    plot(.5:1:3.5,ones(1,4)*.5,'-.k','LineWidth',2);
    ylabel('% correct responses','FontSize',fs)
    xlabel('Phase', 'FontSize',fs)
    title(['Accuracy in ' labs{opt} ], 'FontSize',fs);
    %legend([l1 l2],{'KO','Contr'})
        
    clear Contr KO nContr nKO seContr seKO mtx tmp g
end
set(h1, 'Position', [100, 100, 300, 400]);
keyboard
