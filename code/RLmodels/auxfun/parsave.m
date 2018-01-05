function parsave(fname, model)
    %check if file exists, otherwise change v counter
    [path,name]=fileparts(fname);
    name=[name '.mat'];
    files=dir(path);
    t=structfind(files,'name',name);
    if ~isempty(t)
        num=str2num(fname(end-4));
        if num==9, keyboard; disp('adjust'); end;
        num=num+1;
        fname(end-4)=num2str(num);
    end
    disp([fname '... being saved.'])
    save(fname, 'model');
end 
