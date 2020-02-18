function out = draw_feasibility_dual(model, cone, th, decompose)
    if nargin < 4
        decompose = 1
    end
    N = length(th);
    out.a  = zeros(N, 1);
    out.b  = zeros(N, 1);
    out.x  = cell(N, 1);
    
    if ~isfield(model, 'A')
        A = model.At';
    else
        A = model.A;
    end
    
    if ~isfield(model, 'C')
        C = model.c;
    else
        C = model.C;
    end
    if decompose
        [model_new.A,model_new.b,model_new.C,model_new.K, model_new.info] = ...
            decomposed_subset(A, model.b, C, model.K, cone);
        model_new.pars = model.pars;
    else
        model_new = model;
        model_new.pars.fid = 0;
    end
    %yalmip poses this as a dual optimization problem for some reason
    
    for i = 1:N
        theta = th(i);                
        
        b = zeros(size(model_new.b));

        b(end) = -cos(theta);
        b(end-1) = -sin(theta);                
        %c = c1*cos(theta) + c2*sin(theta);


        %there's a bug here, the chordal decomposition did not add new
        %equality constraints. How to add these in?
        [y,x,info] = sedumi(model_new.A, b, model_new.C, model_new.K, model_new.pars);
        %K = model.K;

        out.a(i) = x(end-1);
        out.b(i) = x(end);
        %x = decomposed_recover(x, model_new.info); 
        %out.x{i} = decomposed_recover(x, model_new.info);
    end
    
    %out.conv = convhull(out.a, out.b);
end

