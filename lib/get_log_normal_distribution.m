function obj = get_log_normal_distribution(cov,M,n) 
V=M*cov;
mu=log(M^2/(sqrt(V+M^2)));
sigma = sqrt(log(V/M^2+1));
% Stiffness matrix for MCS with variation in E
rng(10) ; 
psi = randn(n,1);   %generating normally distributed data
obj.Emc = exp(psi*sigma+mu);   %primary lognormally distributed RV in std normal distribution
