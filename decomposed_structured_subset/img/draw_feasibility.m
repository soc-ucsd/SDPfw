function out = draw_feasibility(model, cone, th, dual)
    N = length(th);
    out.a  = zeros(N, 1);
    out.b  = zeros(N, 1);
    out.x  = cell(N, 1);
    
    if nargin < 4
        dual = 0;
    end
    
    %c = %model.C;
    %c1 = c(1);
    %c2 = c(2);
    
    
%     if dual
%         %testing of DD*
%         [model_new.A,model_new.b,model_new.C,model_new.K, model_new.info] = ...
%             dd_star_convert(model.A, model.b, model.C, model.K);
%         model_new.pars = model.pars;
%     
%     else
        [model_new.A,model_new.b,model_new.C,model_new.K, model_new.info] = ...
            decomposed_subset(model.A, model.b, model.C, model.K, cone);
        model_new.pars = model.pars;
%     end
    
    %yalmip poses this as a dual optimization problem for some reason
    
    for i = 1:N
        theta = th(i);                
        
        if dual
            b = zeros(size(model_new.b));
            
            b(1) = -cos(theta);
            b(2) = -sin(theta); 
            
            [x, y, info] = sedumi(model_new.A, b, model_new.C, model_new.K, model_new.pars);
            %K = model.K;
            %x = decomposed_recover(x, model_new.info);  
            out.a(i) = y(11);
            out.b(i) = y(12);
        else
            C = zeros(size(model_new.C));

            C(1) = -cos(theta);
            C(2) = -sin(theta);                
            %c = c1*cos(theta) + c2*sin(theta);


            %there's a bug here, the chordal decomposition did not add new
            %equality constraints. How to add these in?
            [x,y,info] = sedumi(model_new.A, model_new.b, C, model_new.K, model_new.pars);
            %K = model.K;
            %x = decomposed_recover(x, model_new.info);  
            out.a(i) = x(1);
            out.b(i) = x(2);
        end
%         if length(model.K.s) > 1
%             out.x{i} = {reshape(x(3:11), 3, 3), reshape(x(12:15), 2, 2)};
%         else
%             out.x{i} = reshape(x(3:end), 4, 4);
%         end
    end
    
    %out.conv = convhull(out.a, out.b);
end

