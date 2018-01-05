function logL= getLL(pars,priors,mtx,tp)


%handover parameters
nu=pars(1);     %nu
kappa=pars(2);  %kappa;
b1=pars(3);     %beta1
b2=pars(4);     %beta2
b3=pars(5);     %beta3

days=mtx(:,1);
ndays=length(unique(days));
cond=mtx(:,3);      %Trial types
reward=mtx(:,4);
resp=mtx(:,5);
cue=mtx(:,6);
T=length(cond);
nt1=tp(1);
nt2=tp(2);

%get number of trials per day
realday=unique(days);
for iday=1:ndays
    ind=mtx(:,1)==realday(iday);
    nT(iday)=sum(ind);
end

%initiate qs
qs=zeros(T,4);      %winstay winshift losestay loseshift
qs(1,:)=priors;

%init likelihood
L=zeros(1,T); L(1)=.5;
t=1;

alphat(1:2)=0;
tnewvals=0;

%over all experimental days
for iday=1:ndays
    
    %get phase
    if realday(iday)<nt1, phase=1; b=b1; 
    elseif realday(iday)<nt2 && realday(iday)>=nt1; phase=2; b=b2; 
    elseif realday(iday)>=nt2, phase=3; b=b3; 
    end
    T=nT(iday);
    
    for tt=1:T
        t=t+1;
        if t>sum(nT), continue; end
        if t==tnewvals+1, qs(t-1,:)=qsnew; end
        
        lastrew=reward(t-1);  lastchoice=resp(t-1); lastcue=cue(t-1);
        
        %define responses according to strategies
        stratresp=[lastchoice ~lastchoice lastchoice ~lastchoice]; %winstay winshift losestay loseshift
        
        if lastrew, ind=[1 2]; else ind=[3 4]; end %1winstay 2winshift 3losestay 4loseshift
       
        %get response by strategy
        pr=exp(b*qs(t-1,ind))/sum(exp(b*qs(t-1,ind)));
            
        %get choice and reward at t
        choice=resp(t);
        rew=reward(t);
        
        %which strategy was chosen?
        tmp=choice==stratresp(ind);
        strati=ind(tmp);
        %fill likelihood
        L(t)=pr(tmp);
    
        %delta learning rule
        predR=qs(t-1,strati); delta=rew-predR;
        qs(t,:)=qs(t-1,:);
        qs(t,strati)=qs(t-1,strati)+kappa*alphat(t-1)*delta;

        %update alpha
        alphat(t)=(1-nu)*alphat(t-1) + nu*abs(delta);
        
        %log away: day trial trialperday phase cond reward choice cue q-values
        logdata(t,:)=[realday(iday) t tt phase cond(t) rew choice cue(t) L(t) qs(t,:)];
    end
end
logL=sum(log(L));