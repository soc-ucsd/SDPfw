N = 120;
%correlative sparsity graph of Lehmer-Rosenbrock Function
%Add node to CSP if variables appear together in a monomial
G = eye(N);

%Rosenbrock
%f_R = sum(10*((xi2 + 2*xi1 - xi.^2).^2) + (1-xi -xi3).^2);

G = G + diag(ones(N-1, 1), 1);
G = G + diag(ones(N-2, 1), 2);
G = G + diag(ones(N-3, 1), 3);

G = (G + G')/2;

% for i = 1:N-2
%     
% end



%Lehmer
G(1:N/6, 1:N/6) = 1;



figure(1)
spy(G)

K = cliquesFromSpMatD(G);

lenK = cellfun(@length, K.Set);
ulen = unique(lenK)