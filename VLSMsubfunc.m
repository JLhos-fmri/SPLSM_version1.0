function [t,p,Z,P] = VLSMsubfunc(DAT,GROUP1,GROUP2,COVterm,VALt)
g1 = find(DAT);
group = GROUP2;
group(g1) = GROUP1(g1);
GROUP = term(group);
SLM = 1+GROUP;
contra = [0 1 -1];
for i = 1:size(COVterm,2)
    A = full(COVterm(:,i));
    COVT = term(A);
    SLM = SLM+COVT;
    contra = [contra,0];
end
SLMnew = SurfStatLinMod(VALt,SLM);
SLMres = SurfStatT(SLMnew,contra);
t = SLMres.t;
if t<0
    p = tcdf(SLMres.t,SLMres.df);
else
    p = 1-tcdf(SLMres.t,SLMres.df);
end
[Z, P] = AS_TFRtoZ(t,'T',SLMres.df,[]);
end