clr

load('E:\SiteResp_GMSelect\Data\Dummy motions\Scalar\Scalar.mat');
clearvars -except theta
% theta_mod = theta;
% theta_mod(:,3) = [];

AF = importdata('E:\SiteResp_GMSelect\Data\Dummy motions\Scalar\AmpFac_0.25.txt');
SaR = importdata('E:\SiteResp_GMSelect\Data\Dummy motions\Scalar\SaRock_0.25.txt');

Periods = [0.01, 0.025, 0.05,0.075, 0.1, 0.13, 0.15, 0.18, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6,0.7,0.85,1,1.5,2,2.5,3,3.75,4,5];

No = length(AF); Nt = length(Periods);

for ii = 1:Nt
    
    Pred1(ii,:) = Pred(theta(ii,:), SaR(ii,:));
    
end

res = AF-Pred1;

for ii = 1:Nt
    
   for jj = 1:Nt
       
      corr1(ii,jj) = corr(res(ii,:)',res(jj,:)'); 
       
   end
    
end

corr1 = nearestSPD(corr1);

for ii = 1:Nt
    
   for jj = 1:Nt
       
      COV(ii,jj) = theta(ii,3)*theta(jj,3)*corr1(ii,jj);
       
   end
    
end



% semilogx(Periods,corr1(1,:))
% hold on
% semilogx(Periods,corr1(5,:))
% semilogx(Periods,corr1(10,:))
% semilogx(Periods,corr1(15,:))
% semilogx(Periods,corr1(20,:))
% semilogx(Periods,corr1(24,:))

%% 
% Reshaped to have all three parameters given a period below/above each other
alpha_prior_mu = reshape(theta_mod',Nt*3,1);
alpha_prior_sig = 0.5*eye(length(alpha_prior_mu));
alpha_prop_mu = reshape(theta_mod',Nt*3,1);
alpha_prop_sig = 0.5*eye(length(alpha_prior_mu));

Q = diag(theta(:,3).^2);
nu = 1000;

Nsims = 2000;

alpha_sto(:,1) = mvnrnd(alpha_prop_mu,alpha_prop_sig)';
accept = 0;
for ii = 1:Nsims
    
    alpha_star = mvnrnd(alpha_prop_mu,alpha_prop_sig)';
    res_r = res_vec(reshape(alpha_star',3,Nt)', SaR, Nt, AF);
    Sig_sto(:,:,ii) = iwishrnd((Q+res_r),nu+No);
    
    lnr = Likelihood_vec(reshape(alpha_star',3,Nt)', SaR, Nt, No, AF, reshape(Sig_sto(:,:,ii),Nt,Nt)) + log(mvnpdf(alpha_star, alpha_prior_mu, alpha_prior_sig))...
        -Likelihood_vec(reshape(alpha_sto(:,ii)',3,Nt)', SaR, Nt, No, AF, reshape(Sig_sto(:,:,ii),Nt,Nt)) - log(mvnpdf(alpha_sto(:,ii), alpha_prior_mu, alpha_prior_sig));
    
    if lnr > log(rand)
        alpha_sto(:,ii+1) = alpha_star;
        accept = accept + 1;
    else
        alpha_sto(:,ii+1) = alpha_sto(:,ii);
    end
    
    progressbar(ii/Nsims)
    
    
end
