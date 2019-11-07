function  [Hout, package, time_solve, time_convert] = run_model_star(model, cone, use_mosek)    
    

    %H2 and Hinf norm both require square roots
    if nargin < 3
        use_mosek = 1;
    end    
    
    tic
    
    %Cross over because this is an LMI problem
    %[y,x,infoSeDuMi] = sedumi(-model.A',-model.c,-model.b,model.J,pars);
        
    [A, b, c, K, ~] = decomposed_subset(-model.A',-model.c,-model.b,model.J, cone);
    %pars.fid = 0;
    %pars.fid = 1;    
    %[x, ~, info] = sedumi(A, b, c, K, pars);
    if use_mosek
        prob = sedumi2mosek(A', b, c, K);
        time_convert = toc;
        tic;    
        
        % Set log level (integer parameter)
        %param.MSK_IPAR_LOG = 1;
        % Select interior-point optimizer... (integer parameter)
        param.MSK_DPAR_INTPNT_CO_TOL_REL_GAP = 1.0e-6;
        if all(strcmp(cone, 'dd') )
            param.MSK_IPAR_OPTIMIZER = 'MSK_OPTIMIZER_INTPNT';
            % ... without basis identification (integer parameter)
            param.MSK_IPAR_INTPNT_BASIS = 'MSK_BI_NEVER';            
        end
        
        [r,res] = mosekopt('minimize echo(0)',prob, param);
        %[r,res] = mosekopt('minimize',prob, param)
        
        time_solve = toc;
        if  strcmp(res.sol.itr.prosta, 'PRIMAL_AND_DUAL_FEASIBLE')        
            cost = res.sol.itr.pobjval;
        else
            cost = NaN;
        end

        package = res;
    else
        time_convert = toc;
        tic
        pars.fid = 1;
        %pars.fid = 0;
        [x, y, info] = sedumi(A, b, c, K, pars);
        time_solve = toc;
        cost = c'*x;
        package.x = x;
        package.y = y;
        package.info = info;
    end
    Hout = sqrt(cost);
    %Hout = cost; 
    %Hinf = 
end