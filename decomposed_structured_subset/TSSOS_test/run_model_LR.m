function [cost, package, time_solve, time_convert] = run_model_LR(model, cone, use_mosek, fname)

    if nargin < 3
        use_mosek = 1;
    end           
    
    tic
    
    %Cross over because this is an LMI problem
    %[y,x,infoSeDuMi] = sedumi(-model.A',-model.c,-model.b,model.J,pars);
        
    [A, b, c, K, ~] = decomposed_subset(model.A,model.b,model.c,model.K, cone);

    if use_mosek
        prob = convert_sedumi2mosek(A, b, c, K);
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
        
        [r,res] = mosekopt('minimize',prob, param);
        %if nargin < 4
            %[r,res] = mosekopt('minimize echo(0)',prob, param);
        %else
        %call_str = char(strcat('minimize log(', fname, ')'));
        %[r,res] = mosekopt(call_str,prob, param);
        %end
        %[r,res] = mosekopt('minimize',prob, param);
        
        time_solve = toc;
        if isfield(res, 'sol') && strcmp(res.sol.itr.prosta, 'PRIMAL_AND_DUAL_FEASIBLE')        
            cost = -res.sol.itr.pobjval;
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
        cost = -c'*x;
        package.x = x;
        package.y = y;
        package.info = info;
    end
    %Hout = sqrt(cost);
    %Hout = cost; 
    %Hinf = 
end