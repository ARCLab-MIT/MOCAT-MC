function z=func_dv(Am,s)

% reference: SPACE DEBRIS, Model and Risk Analysis
% eq(3.42)

switch s
    case 'exp'
        mu = 0.2*log10(Am) + 1.85;  % explosion
    case 'col'
        mu = 0.9*log10(Am) + 2.9;   % collision
    otherwise
        warning('exp/col not specified; using explosion parameter')
end


sigma = 0.4;

N = mu+sigma*randn(size(mu));
% N = normrnd(mu,sigma,1,1);

z = 10.^N;        % m/s

