function [model_out, x_fake]= basis_change(x, model, y, dual)
    %Cholesky basis change, based on original paper by Georgina and Amirali
    %
    %Input:
    %   x:      analyitic center of new basis
    %   model:  sedumi format model (struct with fields At, b, c, K)
    %   y:      If using dual variables for outer approximation
    %Output:
    %   model_out: new model under basis change
    
    if nargin < 3
        dual = 0;
        y = 0;
    end
    
    K = model.K;
    
    basis = cell(length(K.s), 1);
    
    Count = K.f + K.l + sum(K.q);
    
    model_out = model;
    x_fake = zeros(size(x));
    
    if dual
        z = model.c - model.At*y;
    end
    
    for i = 1:length(K.s)
        Ksi = K.s(i);
        ind = Count + (1:Ksi^2);
        
        %Find new basis
        x_fake(ind) = reshape(speye(Ksi), [], 1);
        if dual
            Z = reshape(z(ind), Ksi, Ksi);
            Lz = chol(Z + 1e-8*eye(Ksi), 'lower');
            L = inv(Lz);
        else
            X = reshape(x(ind), Ksi, Ksi);        
            L = chol(X, 'lower');
        end
        
        
        basis{i} = L;
        
        %perform basis change
        C = reshape(model.c(ind), Ksi, Ksi);
        C_new = L'*C*L;
        model_out.c(ind) = reshape(C_new, [], 1);
        
        At_curr = model.At(ind, :);
        At_curr0 = At_curr;
        for j = 1:size(model.At, 2)
            At_j = At_curr(:, j);
            At_j_mat = reshape(At_j, Ksi, Ksi);
            
            At_j_mat_new = L'*At_j_mat*L;
            At_curr(:, j) = reshape(At_j_mat_new, [], 1);
        end
        
        model_out.At(ind, :) = At_curr;
        
        Count = Count + Ksi^2;
    end
    
    
    if isfield(model, 'basis')
        for i = 1:length(model.basis)
            basis{i} = model.basis{i}*basis{i};
        end
    else
        model_out.basis = basis;
    end
end