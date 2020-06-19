function [lnAFpred]= Pred(theta, SaR)


a = theta(1);
b = theta(2);
sigma = theta(3);
f3 = theta(4);

f = a+b.*log((SaR+f3)/f3);

lnAFpred = f;

end