function pX = gmmeval(X, mu, sigmasq, wt)
% Evaluate the probability density of data vectors using a GMM model
% Input
%   - X         : N D-dimensional data points (a D-by-N vector)
%   - mu        : mean vectors of Gaussian components (a D-by-K matrix, where
%       each column corresponds to a Gaussian)
%   - sigmasq   : covariance matrices of Gaussian components (a D-by-D-by-K matrix,
%   - wt        : weights of Gaussian components (a 1-by-K vector, sums to 1)
%
% Author: Zhiyao Duan
% Created: 11/1/2012
% Last modified: 11/2/2012

[D, N] = size(X);
K = size(mu,2);

pX = 0;
% sum likelihood for each Gaussian component
for j = 1:K
    if D == 1
        Xzm = X - mu(j);
        pX = pX + wt(j) * (2*pi)^(-1/2) / sqrt(sigmasq(j))...
            * exp(-0.5 * Xzm.*Xzm / sigmasq(j));
    else
        Xzm = X - repmat(mu(:,j), 1, N);
        pX = pX + wt(j) * (2*pi)^(-D/2)/sqrt(det(sigmasq(:,:,j)))...
            * exp(-0.5 * sum(Xzm.*(sigmasq(:,:,j)\Xzm), 1));
    end
end