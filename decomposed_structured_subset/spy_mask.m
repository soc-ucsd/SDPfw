function [mask] = spy_mask(model)
%SPY_MASK generates the sparsity pattern of SDP with At and C
%   Detailed explanation goes here


index = model.K.f + model.K.l + sum(model.K.q);

Cs = model.C(index+1:end);
Ats = model.A(index+1:end, :);

SP = spones(spones(Cs) + sparse(sum(spones(Ats),2)));  % vector of 1s and 0s
mask = reshape(SP, sum(model.K.s), sum(model.K.s));
end

