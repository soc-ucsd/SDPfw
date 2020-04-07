load('Example_SDPfw.mat')

s = svd(full(A));
n_null = nnz(s <= 1e-8);

BlockSize = 4;
pars.fid = 0;

[x0, y0, info0] = sedumi(A, b, c, K, pars);
cost0 = c'*x0;
disp(info0)

% Same process through factorwidth
% opts=struct;
% opts.bfw = 1;
% opts.block = BlockSize;
% 
% [Af, bf, cf, Kf] = factorwidth(A, b, c, K, opts);
% [xf, yf, infof] = sedumi(An, bn, cn, Kn, pars);
% costf = cf'*xf;
%disp(infof)

%Same process through decomposed_subset
[An, bn, cn, Kn, info_dec] = decomposed_subset(A, b, c, K, BlockSize);
[xn, yn, infon] = sedumi(An, bn, cn, Kn, pars);
xnr = decomposed_recover(xn, info_dec);
costn = c'*xnr;
disp(infon)


%eigenvalues of returned solution are in eX
count = 0;
X = cell(length(K.s), 1);
for i = 1:length(K.s)
    Ksi = K.s(i);
    xi = xnr(count + (1:Ksi^2));
    X{i} = reshape(xi, Ksi, Ksi);
    count = count + Ksi^2;
end

eX = cellfun(@eig, X, 'UniformOutput', false);