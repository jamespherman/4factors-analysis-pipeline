function [tbl,chi2stat,pval] = chi2proptest(X)
% [tbl, chi2stat, pval] = chi2proptest(X)
%
% Use the crosstab function to perform a chi-square test on the statistical
% similarity of n proportions a1/A1 vs. a2/A2 vs. ... vs. an/An
%
% The array X should by 2 X n, with "a" values in the first row and "A"
% values in the second row. For example to test 3 proportions:
%
% X = [a1, a2, a3; A1, A2, A3];
%

n = size(X,2);
x1 = [];
x2 = [];

for i = 1:n
    a   = X(1,i);
    A   = X(2,i);
    x1  = [x1; repmat(i,A,1)];
    x2  = [x2; repmat(1,a,1); repmat(2,A-a,1)];
end

[tbl,chi2stat,pval] = crosstab(x1,x2);
end