function [log_likelihood]= Likelihood_vec(theta, SaR, Nt, No, AF, Sigma)

for jj = 1:Nt

    a = theta(jj,1);
    b = theta(jj,2);
    f3 = theta(jj,3);

    f = a+b.*log((SaR(jj,:)+f3)/f3);

    res(jj,:) = log(AF(jj,:))-f;

end

log_likelihood = 0;
for ii = 1:No
    
    log_likelihood = log_likelihood + mvnpdf(res(:,ii),zeros(Nt,1),Sigma);
    
end

end