function H = hessian_weight(VI_dW_dp, weight, N_p, w)
% HESSIAN_WEIGHT - Compute weighted Hessian
%   H = HESSIAN_WEIGHT(VI_DW_DP, WEIGHT, N_P, W)

% Iain Matthews, Ralph Gross
% Carnegie Mellon University, Pittsburgh

if nargin<4 error('Not enough input arguments'); end

H = zeros(N_p, N_p);
for i=1:N_p
  h1 = VI_dW_dp(:,((i-1)*w)+1:((i-1)*w)+w);
  for j=1:N_p
    h2 = VI_dW_dp(:,((j-1)*w)+1:((j-1)*w)+w);
    H(j, i) = sum(sum((weight.*(h1 .* h2))));
  end
end
