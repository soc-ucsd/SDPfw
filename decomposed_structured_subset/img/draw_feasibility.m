function out = draw_feasibility(model, cone, th, DUAL)
    N = length(th);
    out.a  = zeros(N, 1);
    out.b  = zeros(N, 1);
    out.x  = cell(N, 1);
    out.y  = cell(N, 1);
    
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
    
    [model_new.A,model_new.b,model_new.C,model_new.K, model_new.info] = ...
        decomposed_subset(A, model.b, C, model.K, cone, DUAL);
    model_new.pars = model.pars;
 
    %yalmip poses this as a dual optimization problem for some reason
    
    for i = 1:N
        theta = th(i);                
        
        C = zeros(size(model_new.C));

        C(1) = -cos(theta);
        C(2) = -sin(theta);                
        %c = c1*cos(theta) + c2*sin(theta);


        %there's a bug here, the chordal decomposition did not add new
        %equality constraints. How to add these in?
        [x,y,info] = sedumi(model_new.A, model_new.b, C, model_new.K, model_new.pars);
        %K = model.K;

        out.a(i) = x(1);
        out.b(i) = x(2);
        %x = decomposed_recover(x, model_new.info); 
        out.x{i} = decomposed_recover(x, model_new.info);
        out.y{i} = y(1:length(model.b));
    end
    
    %out.conv = convhull(out.a, out.b);
end

