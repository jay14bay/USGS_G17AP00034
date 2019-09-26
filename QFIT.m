function f=QFIT(beta,X)

b1=beta(1);
b2=beta(2);
freq=X(:,1);

f=log(b1) + b2.*log(freq);

end