function [res_r]= res_vec(theta, SaR, Nt, AF)

for jj = 1:Nt

    a = theta(jj,1);
    b = theta(jj,2);
    f3 = theta(jj,3);

    f = a+b.*log((SaR(jj,:)+f3)/f3);

    res(jj,:) = log(AF(jj,:))-f;

end

res_r = res*res';

end