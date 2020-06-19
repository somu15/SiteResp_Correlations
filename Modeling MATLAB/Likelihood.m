function [L,lnAFpred]= Likelihood(theta, SaR, AF)


a = theta(1);
b = theta(2);
sigma = theta(3);
f3 = theta(4);

f = a+b.*log((SaR+f3)/f3);

lnAFpred = f;

D=(log(AF)-lnAFpred);

x = log(1/sigma*exp(-.5*((D)/sigma).^2));

L = -sum(x);


end