function [ztau ptau] = kendalltau(N,tau)
% 
% Given the number of data points N for a given number of spins and transfer
% and given the Kendall tau, this function computes the statistic Ztau and
% the confidence p
% 
% Standard deviation of Kendall tau
%
stdv=sqrt(2.*(2.*N+5)./(9.*N.*(N-1)))
%
% Kendall statistics
%
ztau=tau./stdv;
%
% Compute p from normal distribution
%
ptau=1-0.5*erfc(-ztau/realsqrt(2))
%
end

