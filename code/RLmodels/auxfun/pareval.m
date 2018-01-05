function ratname = pareval(ratname, x)
    eval(sprintf([ratname '=struct([])']));
    eval([ratname '=x']);
    ratname=eval(ratname);
end