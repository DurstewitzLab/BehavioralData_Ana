function make_figs_4A()
%plots strategy frequencies

cd ..
pati=[pwd '/data/']
addpath([pwd '/code/helpfunctions/']);

dum1=strcat(pati,'Group1*.mat');    %Knockout group (Cav1.2NesCre)
dum2=strcat(pati,'Group2*.mat');    %Control (Cav1.2flfl)
files1=dir(dum1); files2=dir(dum2);
files=[files1; files2];
[nFiles, ~]=size(files);


%loop over files 
for i=1:nFiles
    file=[pati files(i).name];
    disp(file)
    f=files(i).name; 
    vpn=f(end-5:end-4);
    vpn(vpn=='_')=[];
    vpn=str2num(vpn);
    grouplabs={'group1','group2'}; %KO Contr
    
    dat=load([pati files(i).name]);  
    mtx=dat.mtx;
    
    %----------------------------------------------------------------------
    res(i).rat=vpn;
    res(i).group=dat.g;
    res(i).grouplab= grouplabs{dat.g};
    res(i).phasestart=dat.tp;
    res(i).mtx=dat.mtx;

    %shown win-shift/stay and lose-shift/stay strategies
    eval(['res(i).p_winstay=dat.stratvals(:,1);'])
    eval(['res(i).p_winshift=dat.stratvals(:,2);'])
    eval(['res(i).p_losestay=dat.stratvals(:,3);'])
    eval(['res(i).p_loseshift=dat.stratvals(:,4);'])
    
end
clearvars -except pat* res

%create figure on frequency of strategies
createfig(res);


%--------------------------------------------------------------------------
function createfig(data)

xx=[1 109 124 143];
n=size(data,2);
labs={'p_winstay','p_winshift','p_losestay','p_loseshift'};

fs=16;
h1= figure('color','white'); hold on; box on
h2= figure('color','white'); hold on;  box on
c={'r','m','b','c'};

%sort the data
strat=[2 4];
for opti=strat%length(labs)
    opt=opti;
    mtx=nan(24,143);
    
    m=0;
    for i=1:n
        
        tmp=data(i); g(i)=tmp.group;
        y1=1:tmp.phasestart(1)-1;
        y2=tmp.phasestart(1):tmp.phasestart(2)-1;
        y3=tmp.phasestart(2):tmp.phasestart(3);
        
        eval(['dum1=tmp.' labs{opt} '(y1);'])
        eval(['dum2=tmp.' labs{opt} '(y2);'])
        eval(['dum3=tmp.' labs{opt} '(y3);'])
        
        x1=y1; x2=xx(2):xx(2)+length(y2)-1; x3=xx(3):xx(3)+length(y3)-1;
        mtx(i,x1)=dum1;
        mtx(i,x2)=dum2;
        mtx(i,x3)=dum3;
        
        clear dum* y*
    end

    %create the figure
    cols=[.98 .98 .98; .95 .95 .95;  .92 .92 .92];
    Contr=mtx(g==2,:); nContr=sum(~isnan(Contr)); mContr=nanmean(Contr); seContr=nanstd(Contr)./sqrt(nContr);
    KO=mtx(g==1,:); nKO=sum(~isnan(KO)); mKO=nanmean(KO); seKO=nanstd(KO)./sqrt(nKO);
            
    
    h=colormap('jet'); c1=h(5,:); c2=h(end-6,:);
    figure(h1); hold on
    for iphase=1:3
        x=xx(iphase):xx(iphase+1)-1;
        if opti==strat(1)
            h=area([x(1) x(end)+1],[.995 .995],'FaceColor',cols(iphase,:),'EdgeColor','none');
            hpatch = get(h,'children'); set(hpatch,'FaceAlpha',0.8)
        end
        
        xin=~isnan(mContr(x));
        if sum(xin)>0
            l1(opti)=shadedErrorBar(x(xin),mContr(x(xin)),seContr(x(xin)),{'color',c{opti},'LineWidth',2})
        end
    end
    set(gca,'FontSize',fs); xlim([0 143])
    plot(xx(1):xx(end),.5,'-.k','LineWidth',2); str=labs{opt};
    text(5,.9, 'Ca_v1.2^{fl/fl}','FontSize',20,'FontWeight','bold')
    ylabel('relative frequency','FontSize',fs);
    xlabel('day #', 'FontSize',fs)    
    
    %Knockout group
    figure(h2);hold on
    h=colormap('jet'); c1=h(5,:); c2=h(end-6,:);
    for iphase=1:3
        x=xx(iphase):xx(iphase+1)-1;
        if opti==strat(1) 
            h=area([x(1) x(end)+1],[.995 .995],'FaceColor',cols(iphase,:),'EdgeColor','none');
            hpatch = get(h,'children'); set(hpatch,'FaceAlpha',0.8)
        end
        
        if sum(xin)>0
            xin=~isnan(mKO(x));
            l2(opti)=shadedErrorBar(x(xin),mKO(x(xin)),seKO(x(xin)),{'color',c{opti},'LineWidth',2})
        end        
    end
    set(gca,'FontSize',fs); xlim([0 143])
    plot(xx(1):xx(end),.5,'-.k','LineWidth',2); str=labs{opt};
    text(5,.9, 'Ca_v1.2^{NesCre}','FontSize',20,'FontWeight','bold')
    ylabel('relative frequency','FontSize',fs);
    xlabel('day #', 'FontSize',fs)
    clear Contr KO nContr nKO seContr seKO mtx tmp g
end
figure(h1); legend([l1.mainLine,l1.mainLine],'winshift','loseshift')
figure(h2); legend([l2.mainLine,l2.mainLine],'winshift','loseshift')

    keyboard

