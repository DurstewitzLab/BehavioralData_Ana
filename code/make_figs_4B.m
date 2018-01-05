function make_figs_4B()
%collects and plots reward probabilities for different strategies

clear; close all; clc

cd ..
pati=[pwd '/data/']
addpath([pwd '/code/helpfunctions/']);
%uses function shadedErrorBar written by Rob Campbell - November 2009

dum1=strcat(pati,'Group1*.mat');    %Knockout group (Cav1.2NesCre)
dum2=strcat(pati,'Group2*.mat');    %Control (Cav1.2flfl)
files1=dir(dum1); files2=dir(dum2);
files=[files1; files2];
[nFiles, ~]=size(files);

xx=[1 109 124];
for i=1:nFiles
    file=[pati files(i).name];
    disp(file)
    f=files(i).name;
    vpn=f(end-5:end-4);
    vpn(vpn=='_')=[];
    vpn=str2num(vpn);
    
    dat=load(file);
    g=dat.g;
    tp=dat.tp;
    group(i)=g; %2=Contr, 1=KO
    data=dat.mtx;
    ndays=length(unique(data(:,1)));
    
    nstrats=zeros(ndays,4); nrewards=zeros(ndays,4);
    for iday=1:ndays
        ind=data(:,1)==iday;
        datad=data(ind,:);
        T=size(datad,1);
        for t=2:T
            win=datad(t-1,4); rew=datad(t,4); shiftresp=datad(t,5)-datad(t-1,5);
            switch win
                case 0
                    if shiftresp==0 %losestay
                        index=3;
                    else            %loseshift
                        index=4;
                    end
                case 1
                    if shiftresp==0 %winstay
                        index=1;
                    else            %winshift
                        index=2;
                    end
            end
            nstrats(iday,index)=nstrats(iday,index)+1;
            nrewards(iday,index)=nrewards(iday,index)+rew;
        end
    end
    t1=1:tp(1)-1; t2=tp(1):tp(2)-1; t3=tp(2):tp(3);
    
    %sort into overall matrix per phase
    pmtxwinstay(i,:)=getV(tp,xx,nrewards(:,1)./nstrats(:,1));
    pmtxwinshift(i,:)=getV(tp,xx,nrewards(:,2)./nstrats(:,2));
    pmtxlosestay(i,:)=getV(tp,xx,nrewards(:,3)./nstrats(:,3));
    pmtxloseshift(i,:)=getV(tp,xx,nrewards(:,4)./nstrats(:,4));
    
    clear nstrats nrewards  dat index  tp
end

createfig(pmtxwinstay,pmtxwinshift,pmtxlosestay,pmtxloseshift,group)



% %--------------------------------------------------------------------------
function outV=getV(tp,xx,inV)

t1=1:tp(1)-1; t2=tp(1):tp(2)-1; t3=tp(2):tp(3);
tmp=nan(1,143);
tmp(t1)=inV(t1);
tmp(xx(2):xx(2)+length(t2)-1)=inV(t2);
tmp(xx(3):xx(3)+length(t3)-1)=inV(t3);
outV=tmp;



%--------------------------------------------------------------------------
function createfig(pmtxwinstay,pmtxwinshift,pmtxlosestay,pmtxloseshift,g)
xx=[1 109 124 143];
n=size(pmtxwinstay,2);
labs={'pmtxwinstay','pmtxwinshift','pmtxlosestay','pmtxloseshift'};
fs=16;
h1= figure('color','white'); hold on; box on
h2= figure('color','white'); hold on;  box on
c={'r','m','b','c'};

%sort the data
strat=[1 2 3 4];
for opti=strat
    
    eval(['mtx=' labs{opti} ';'])
    
    %create the figure
    cols=[.98 .98 .98; .95 .95 .95;  .92 .92 .92];
    Contr=mtx(g==2,:); nContr=sum(~isnan(Contr)); mContr=nanmean(Contr); seContr=nanstd(Contr)./sqrt(nContr);
    KO=mtx(g==1,:); nKO=sum(~isnan(KO)); mKO=nanmean(KO); seKO=nanstd(KO)./sqrt(nKO);
    
    
    %Contr
    h=colormap('jet'); c1=h(5,:); c2=h(end-6,:);
    figure(h1); hold on
    for iphase=1:3
        x=xx(iphase):xx(iphase+1)-1;
        if opti==strat(1)
            h=area([x(1) x(end)+1],[1 1],'FaceColor',cols(iphase,:),'EdgeColor','none');
            hpatch = get(h,'children'); set(hpatch,'FaceAlpha',0.8);
        end
        
        xin=~isnan(mContr(x));
        if sum(xin)>0
            l1(opti)=shadedErrorBar(x(xin),mContr(x(xin)),seContr(x(xin)),{'color',c{opti},'LineWidth',2});
        end
    end
    set(gca,'FontSize',fs); xlim([0 143])
    plot(xx(1):xx(end),.5,'-.k','LineWidth',2); str=labs{opti};
    text(5,.9, 'Ca_v1.2^{fl/fl}','FontSize',20,'FontWeight','bold');
    ylabel('reward probability','FontSize',fs);
    xlabel('day #', 'FontSize',fs)
    ylim([0 1.005])
    
    %KO
    figure(h2);hold on
    h=colormap('jet'); c1=h(5,:); c2=h(end-6,:);
    for iphase=1:3
        x=xx(iphase):xx(iphase+1)-1;
        if opti==strat(1)
            
            h=area([x(1) x(end)+1],[1 1],'FaceColor',cols(iphase,:),'EdgeColor','none');
            hpatch = get(h,'children'); set(hpatch,'FaceAlpha',0.8)
        end
        
        if sum(xin)>0
            xin=~isnan(mKO(x));
            l2(opti)=shadedErrorBar(x(xin),mKO(x(xin)),seKO(x(xin)),{'color',c{opti},'LineWidth',2});
        end
    end
    set(gca,'FontSize',fs); xlim([0 143])
    plot(xx(1):xx(end),.5,'-.k','LineWidth',2); str=labs{opti};
    text(5,.9, 'Ca_v1.2^{NesCre}','FontSize',20,'FontWeight','bold')
    ylim([0 1.005])
    
    ylabel('reward probability','FontSize',fs);
    xlabel('day #', 'FontSize',fs)
    clear Contr KO nContr nKO seContr seKO mtx tmp
end
figure(h2);
legend([l2.mainLine,l2.mainLine],{'winstay','winshift','losestay','loseshift'},'FontSize',14);
keyboard
