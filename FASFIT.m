function f=FASFIT(beta,X)

b1=beta(1);
b2=beta(2);
R=X(:,1);

f=b1 + b2.*(R) - 0.5.*log(R);
end