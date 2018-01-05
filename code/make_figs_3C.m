function make_figs_2C()

clear; close all; clc

cd ..
pati=[pwd '/data/']

dum1=strcat(pati,'Group1*.mat');    %Knockout group (Cav1.2NesCre)
dum2=strcat(pati,'Group2*.mat');    %Control (Cav1.2flfl)
files1=dir(dum1); files2=dir(dum2);
files=[files1; files2];
[nFiles, ~]=size(files);


grouplabs={'Group1','Group2'}; %KO vs. Control

for i=1:nFiles
    file=[pati files(i).name];
    disp(file)
    f=files(i).name;
    vpn=f(end-5:end-4);
    vpn(vpn=='_')=[];
    vpn=str2num(vpn);
    
    %load file, and check for errors
    dat=load([pati files(i).name]);
    
    mtx=dat.mtx;
    days=unique(mtx(:,1));      %no days
    %----------------------------------------------------------------------
    %get information on correct trials etc
    res(i).rat=vpn;
    res(i).group=dat.g;
    res(i).grouplab= grouplabs{dat.g};
    res(i).phasestart=dat.tp;
    res(i).mtx=dat.mtx;
    tp=dat.tp;
    
    %%%average data
    ShiftLC2=[]; ShiftHC2=[]; ShiftLC3=[]; ShiftHC3=[]; StayLC3=[]; StayHC3=[];
    for iday=tp(1):tp(3)
        ind=mtx(:,1)==iday;
        mtxd=mtx(ind,:);
        for icond=1:5
            ind=mtxd(:,3)==icond;
            switch icond
                case 1   %only in third phase
                    if iday>=tp(2)
                        StayHC3=[StayHC3 sum(mtxd(ind,4))/sum(ind)];
                    end
                case 2
                    if iday>=tp(2)
                        StayLC3=[StayLC3 sum(mtxd(ind,4))/sum(ind)];
                    end
                case 3
                    if iday>=tp(2)
                        ShiftHC3=[ShiftHC3 sum(mtxd(ind,4))/sum(ind)];
                    else
                        ShiftHC2=[ShiftHC2 sum(mtxd(ind,4))/sum(ind)];
                    end
                case 4 %ambig trial..
                    
                case 5
                    if iday>=tp(2)
                        ShiftLC3=[ShiftLC3 sum(mtxd(ind,4))/sum(ind)];
                    else
                        ShiftLC2=[ShiftLC2 sum(mtxd(ind,4))/sum(ind)];
                    end
            end
        end
    end
    
    results(i,:)=[vpn dat.g mean(ShiftHC2) mean(ShiftHC3) mean(StayHC3) mean(ShiftLC2) mean(ShiftLC3) mean(StayLC3)];
    
    clearvars -except dat nFiles files pat* res* d i* mtx vpn g grouplabs
end

clearvars -except pat* res* tresults
createfig(results);




%--------------------------------------------------------------------------
function createfig(data)

Contr=data(data(:,2)==2,:); KO=data(data(:,2)==1,:); KO(KO(:,1)==17,:)=[];
h=colormap('jet'); c1=h(5,:); c2=h(end-6,:); %c1=Controls c2=Knockouts
close

%dependent t-tests for phase II 
y1_P2HC = KO(:,3); 
y1_P2LC=KO(:,6); 
y2_P2HC=Contr(:,3);
y2_P2LC=Contr(:,6);



% [~, p1]=ttest(y1_P2HC,y1_P2LC); 
% [~, p2]=ttest(y2_P2HC,y2_P2LC);

%dependent t-tests for phase III
y1_P3HC=mean(KO(:,4:5),2); 
y2_P3HC=mean(Contr(:,4:5),2);
y1_P3LC=mean(KO(:,7:8),2); 
y2_P3LC=mean(Contr(:,7:8),2);



% [~, p1]=ttest(y1_P3HC,y1_P3LC); 
% [~, p2]=ttest(y2_P3HC,y2_P3LC);

fs=8;
h1=figure('color','white'); hold on
errorbar(1,nanmean(y1_P2HC),nanstd(y1_P2HC)/sqrt(length(y1_P2HC)),'-k');
l1=bar(1,nanmean(y1_P3HC),'FaceColor','k','EdgeColor','k');
bar(1,nanmean(y1_P2HC),'FaceColor',c2, 'EdgeColor',c2);
errorbar(2,nanmean(y1_P2LC),nanstd(y1_P2LC)/sqrt(length(y1_P2LC)),'-k');
l2=bar(2,nanmean(y1_P3LC),'FaceColor','w', 'EdgeColor','k')
bar(2,nanmean(y1_P2LC),'FaceColor','w', 'EdgeColor',c2);
errorbar(4,nanmean(y2_P2HC),nanstd(y2_P2HC)/sqrt(length(y2_P2HC)),'-k');
bar(4,nanmean(y2_P2HC),'FaceColor',c1, 'EdgeColor',c1);
errorbar(5,nanmean(y2_P2LC),nanstd(y1_P2LC)/sqrt(length(y2_P2LC)),'-k');
bar(5,nanmean(y2_P2LC),'FaceColor','w', 'EdgeColor',c1);
xlim([0 6]); ylim([.5 .9]);
set(gca,'XTick',[1.5 4.5],'XTickLabel','', 'FontSize',fs, 'YTick',[0 :.1: 1],'YTickLabel',[0:10:100]);
plot(0:6, ones(1,7)*.5,'k');
ylabel('% correct responses','FontSize',fs)
xlabel('Phase II','FontSize',fs); box on
legend([l1 l2],{'High contrast','Low contrast'},'Location','Best')
set(h1, 'Position', [100, 100 180 180]);




%phase III
h1=figure('color','white'); hold on
errorbar(1,nanmean(y1_P3HC),nanstd(y1_P3HC)/sqrt(length(y1_P3HC)),'-k');
bar(1,nanmean(y1_P3HC),'FaceColor',c2,'EdgeColor',c2);
errorbar(2,nanmean(y1_P3LC),nanstd(y1_P3LC)/sqrt(length(y1_P3LC)),'-k');
bar(2,nanmean(y1_P3LC),'FaceColor','w', 'EdgeColor',c2);
errorbar(4,nanmean(y2_P3HC),nanstd(y2_P3HC)/sqrt(length(y2_P3HC)),'-k');
bar(4,nanmean(y2_P3HC),'FaceColor',c1, 'EdgeColor',c1);
errorbar(5,nanmean(y2_P3LC),nanstd(y1_P3LC)/sqrt(length(y2_P3LC)),'-k');
bar(5,nanmean(y2_P3LC),'FaceColor','w', 'EdgeColor',c1);
xlim([0 6]); ylim([.5 .9]);
set(gca,'XTick',[1.5 4.5],'XTickLabel','', 'FontSize',fs, 'YTick',[0 :.1: 1],'YTickLabel',[0:10:100]);
plot(0:6, ones(1,7)*.5,'k');
xlabel('Phase III','FontSize',fs);
ylabel('% correct responses','FontSize',fs);
set(h1, 'Position', [100, 100 180 180]); box on


keyboard

