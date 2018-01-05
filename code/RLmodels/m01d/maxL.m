function maxL(mname,par0, lb,ub,pati)
%maximize log likelihood


[~, modelname]=fileparts(mname(1:end-1));
mno=str2num(modelname(3));

version='_v01';
pato=[mname 'MLE/'];
if ~isdir(pato),    mkdir(pato); end

dum1=strcat(pati,'Group1*.mat');    %Knockout group (Cav1.2NesCre)
dum2=strcat(pati,'Group2*.mat');    %Control (Cav1.2flfl)
files1=dir(dum1); files2=dir(dum2);
files=[files1; files2];
[nFiles, ~]=size(files);

options = optimset('Display','off','MaxFunEvals',500,...
    'TolFun',1e-06,...
    'Algorithm','interior-point');

parInit=par0;
grouplabs={'group1','group2'};

%pooli=parpool('local',12); %for parallel computing
%par
for i=1:nFiles
    filename=files(i).name;
    dat=load([pati filename]);
    
    strats=dat.stratvals(1,:);  %strategy frequency
    tp=dat.tp;                  %how many days should simulation run
    data=dat.mtx(:,1:6);        
    sample=grouplabs{dat.g};
    
    vpn=filename(end-5:end-4);
    vpn(vpn=='_')=[];
    vpn=str2num(vpn);
    
    
    Mtx=ones(1,10)*-10000;
    zsp=zeros(1,3+size(parInit,2)*2);
    A=[];B=[];Aeq=[]; Beq=[];
    for j = 1:size(parInit,1) %loop through init conditions
        par0=parInit(j,:);
        p=fmincon(@(p) -getLL(p,strats,data,tp),par0,A,B,Aeq,Beq,lb,ub); %maximize logL
        l=getLL(p,strats,data,tp)
        disp(['Estimated par:' num2str(p) ' LL:' num2str(l)])
        tmp=[vpn dat.g par0 p l];
        if  l>Mtx(end), Mtx=tmp; end
        zsp=[zsp; tmp];
    end
    [~,ind]=sort(zsp(:,end),'descend'); zsp=zsp(ind,:);
        
    %save all info to mat-file
    vpname=[modelname '_' num2str(vpn) version];
    npars=size(parInit,2);
    x=struct('model',vpname);
    x.model_no=mno;
    x.winpar=Mtx(1,(npars+3):(2*npars+2));
    x.winparnames={'nu','kappa','beta1','beta2','beta3'};
    x.LL=Mtx(1,end);
    x.npar=npars;
    x.allestimates=zsp;
    x.estimationoptions=options;
    x.ub=ub;
    x.lb=lb;
    x.pathdata=pati;
    x.pathmodel=pwd;
    x.pathoutput=pato;
    x.info.vp=vpn;
    x.info.dat=dat;
    x.info.sample=sample;
    x.info.group=dat.g;
    varname=pareval(vpname,x);
    parsave([pato vpname '.mat'], varname);
end
delete(pooli)
