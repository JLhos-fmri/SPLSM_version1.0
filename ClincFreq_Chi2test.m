function [p_ANOVA,Q_ANOVA,p_two,q_two,or_two,or_twolow,or_twohigh] = ClincFreq_Chi2test(x,C)
[p_ANOVA,Q_ANOVA] = chi2test(x);
for i = 1:size(C,1)
    X = [x(C(i,1),:);x(C(i,2),:)];
    [p_two(i,1),q_two(i,1)] = chi2test(X);
    freqv = x(C(i,1),:)./x(C(i,2),:);
    if freqv(1)<freqv(2)
        q_two(i,1) = q_two(i,1)*-1;
    end
    [or_two(i,1),or_twolow(i,1),or_twohigh(i,1)] = Chi2testOR(X(1),X(2),X(3),X(4));
end

end