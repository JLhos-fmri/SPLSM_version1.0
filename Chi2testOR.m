function [orval orcilow orcihigh]=Chi2testOR(A,B,C,D)
if B*C~=0
    orval = (A*D)/(B*C);
    lnorcilow = log(orval)-1.96*sqrt(1/A+1/B+1/C+1/D);
    lnorcihigh = log(orval)+1.96*sqrt(1/A+1/B+1/C+1/D);
    orcilow = exp(lnorcilow);
    orcihigh = exp(lnorcihigh);
else
    orval = NaN;
    orcilow = NaN;
    orcihigh = NaN;
end

end