function [accTot,accShift,accStay]=sim_cue(stratvals1,stratvals2,tp,nT)

close all; clc

nt1=tp(1);
nt2=tp(2);
ndays=size(stratvals1,1);

cond(1:2)=[1 1];
c(1:2)=[1 1 ];      %assume first response correct
resp(1:2)=[1 1];    %0=left 1=right (right was correct on first condition 1)
cue(1:2)=[1 1];

t=2;
for iday=1:ndays
    if iday<nt1, task=1;
    elseif iday<nt2 && iday>=nt1; task=2;
    elseif iday>=nt2, task=3;
    end
    if nt2<nt1 && iday>=nt1, break; end
    
    counton=0;
    countsw=0;
    T=nT(iday);
    for tt=1:T
        t=t+1;
        
        lastrew=c(t-1); lastcue=cue(t-1);
        [cond(t), cue(t)]=getCond(task,lastrew,counton,cond,t,lastcue);
        if cond(t)==1 || cond(t)==3 %switch between probs for low/high contrast
            stratvals=stratvals1;
        else
            stratvals=stratvals2;
        end
        
        if cue(t)==1, corrresp=1; else corrresp=0; end
                
        if cue(t)==1, ind=[1 2]; noind=[3 4]; else ind=[3 4]; noind=[1 2]; end %cue1resp1 cue1resp0 cue2resp1 cue2resp0
        
        %get probability and draw strategy
        choiceprob=stratvals(iday,ind);
        if sum(isnan(choiceprob))>0,
            choiceprob=stratvals(iday,noind);
        end
        pr=choiceprob;
        strat(t)=randsample(ind,1,true,pr);

        %get response according to strategy
        switch strat(t)
            case 1 %cue1r1  (always correct)
                resp(t)=1;
            case 2 %cue1r0
                resp(t)=0;
            case 3 %cue2r1  
                resp(t)=1;  
            case 4 %cue2r0  (always correct)
                resp(t)=0;
        end
        
        %reward yes/no
        if resp(t)==corrresp, rew=1; c(t)=1; else rew=0; c(t)=0; end
        if cond(t)==4, rew=1; c(t)=1; end %ambiguous always correct
                
        %keep track of trial numbers
        if cond(t)==1 && c(t)==1, counton=counton+1; end
        if cond(t)~=1, counton=0; end
        if cond(t)~=3, countsw=0; end
        if cond(t)==3 && c(t)==1, countsw=countsw+1; end
        
        log(t,:)=[iday t tt cond(t) c(t) resp(t) cue(t) task];
    end
    
    %accuracies for different trial types
    ind=log(:,1)==iday;
    accTot(iday)=sum(c(ind)==1)/sum(ind);
    tmp=log(ind,:);
    ind=tmp(:,4)==1 | tmp(:,4)==2;
    accStay(iday)=sum(tmp(ind,5))/sum(ind);
    ind=tmp(:,4)==3 | tmp(:,4)==5;
    accShift(iday)=sum(tmp(ind,5))/sum(ind);
    
    ind=tmp(:,4)==1 ;
    accStayHC(iday)=sum(tmp(ind,5))/sum(ind);
    ind=tmp(:,4)==2 ;
    accStayLC(iday)=sum(tmp(ind,5))/sum(ind);
    ind=tmp(:,4)==3 ;
    accShiftHC(iday)=sum(tmp(ind,5))/sum(ind);
    ind=tmp(:,4)==5 ;
    accShiftLC(iday)=sum(tmp(ind,5))/sum(ind);
end


%--------------------------------------------------------------------------
%draw next condition and cue position for simulation
function [nextcond, nextcue]=getCond(task,lastrew,counton,cond,t,lastcue)

%probabilities per task
p1=[.5 0 .5 0 0];
p2=[.5 0 .25 0 .25];
p3=[.5 .125 .125 .125 .125];
eval(['p=p' num2str(task) ';']);

condlast=cond(t-1);
if t>2
    condlastlast=cond(t-2);
else
    condlastlast=10;
end

newcue=1:2;
newcue(newcue==lastcue)=[];


switch task
    case 1 %phase 1 (50/50) and correction
        if ~lastrew   %correction trial
            nextcond=condlast;
            nextcue=lastcue;
        elseif lastrew && counton>=3
            nextcond=3;
            nextcue=newcue;
        else
            nextcond=randsample(1:5,1,true,p);
            if nextcond==1, nextcue=lastcue; else nextcue=newcue; end
        end
        
    case 2 %phase 2
        if ~lastrew
            nextcond=condlast;
            nextcue=lastcue;
        end
        
        
        if lastrew && counton < 5
            if condlast~=1
                nextcond=1;
                nextcue=lastcue;
            elseif lastrew && condlast==1 && condlastlast~=1 %minimum of 2 ongoings
                nextcond=1;
                nextcue=lastcue;
            else
                nextcond=randsample(1:5,1,true,p);
                if nextcond==1 || nextcond==2|| nextcond==4 , nextcue=lastcue; else nextcue=newcue; end
            end
        elseif lastrew && counton>=5
            nextcond=randsample([3 5],1,true,[.5 .5]); %dann wird cue gewechselt
            nextcue=newcue;
        end
        
    case 3  %phase 3
        if counton>=3, %maximum number of ongoing trials
            nextcond=randsample(2:5,1,true,[.25 .25 .25 .25]);
            if nextcond==1||nextcond==2, nextcue=lastcue;
            elseif nextcond==4,
                nextcue=randsample([1 2], true, 1, [.5 .5]);
            else
                nextcue=newcue;
            end
        end
        
        if counton~=3
            if condlast~=1
                nextcond=1; nextcue=lastcue;
            elseif sum(cond(t-3:t-1))==3  
                nextcond=randsample(2:5,1,true,[.25 .25 .25 .25]);
                if nextcond==1||nextcond==2, nextcue=lastcue;
                elseif nextcond==4,
                    nextcue=randsample([1 2], true, 1, [.5 .5]);
                else
                    nextcue=newcue;
                end
                
            else
                nextcond=randsample(1:5,1,true,p);
                if nextcond==1||nextcond==2, nextcue=lastcue;
                elseif nextcond==4,
                    nextcue=randsample([1 2], true, 1, [.5 .5]);
                else
                    nextcue=newcue;
                end
            end
        end
end
