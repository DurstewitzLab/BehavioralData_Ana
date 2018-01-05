%estimate RL models
clear all; close all; clc

addpath([pwd '/auxfun/']);

optionall=false; %all files or specific files
model={'m02d','m02d'};
datapath=[fileparts(fileparts(pwd)) '/data/'];


%get all subdirectories for model
if optionall
    files = dir;
    filenames = {files.name};
    subdirs = filenames([files.isdir]);
    zsp=strmatch('fs_m',subdirs);
    for s = 1:length(zsp)
        subdir{s,:} = [pwd '/' subdirs{zsp(s)}];
    end
    clear filenames files s subdirs
else
    for s=1:length(model)
        subdir{s,:}=[pwd '/' model{s}];
    end
end

%%%-> boundary conditions (alpha, nu, beta1, beta2, beta3)
malllb=[0 0 0 0 0 ];
mallub=[1 1 500 500 500];
for ibound=1:size(mallub,1)     %for all boundary conditions
    alllb=malllb(ibound,:);
    allub=mallub(ibound,:);
    
    for idir=1:size(subdir,1)   %loop over all model folders
        npwd=char(subdir(idir,:));
        cd(npwd);
        modelpath=[npwd '/'];

        %determine parameters of model
        parinfofile=[npwd '/parameters.txt'];
        fid=fopen(parinfofile);
        c=textscan(fid,'%s %s %s %s %s');
        for ipar=1:size(c,2)
            dum=c{ipar};
            parlabels{ipar}=dum(1);
            parsused(ipar)=str2num(dum{2});
        end
        fclose(fid);
        parsused=parsused>0;
        
        %set initial conditions 
        allpar0 = {[0:.1:1],[0:.1:1],[0:.2:3],[0:.2:3],[0:.2:3]};
        
        modelpar0=allpar0(parsused);
        lb=alllb(parsused);
        ub=allub(parsused);
        npar=size(modelpar0,2);
        switch npar
            case 3
                [x,y,z]=ndgrid(modelpar0{:}); par0=[x(:),y(:),z(:)];
            case 4
                [x,y,z,w]=ndgrid(modelpar0{:}); par0=[x(:),y(:),z(:),w(:)];
            case 5
                [x,y,z,w,l]=ndgrid(modelpar0{:}); par0=[x(:),y(:),z(:),w(:),l(:)];
            otherwise
                display('wrong number of parameters, please check'); keyboard;
        end
        

        %maximize log likelihood and save info to file
        maxL(modelpath,par0,lb,ub,datapath);  
       
        %protocol
        s=datestr(now); protpath=fileparts(modelpath(1:end-1));
        txt=sprintf('MODEL ESTIMATION DONE in %s  at: %s\n',modelpath,s);
        protocol=[protpath '/protocol.txt'];
        fid=fopen(protocol,'at'); fprintf(fid,'%s',txt); fclose(fid);

        
        clearvars -except subdir idir all* mall* iinit*
    end
end