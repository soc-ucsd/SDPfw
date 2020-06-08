function [Hout, time_solve, sdp_opt] = sea_star_reference(model, QUIET, use_mosek)
    A = model.A;
    b = model.b; 
    c = model.c;
    K = model.K;
    
    if nargin < 2
        QUIET = 1;
    end
    
    if nargin < 3
        use_mosek = 1;
    end


    if use_mosek
        %prob = convert_sedumi2mosek(A', b, c, K);
        %prob0 = convert_sedumi2mosek(-model.A',-model.c,-model.b,model.J);
        %[r0, res0] = mosekopt('minimize', prob0);       
        prob = sedumi2mosek(A', b, c, K);
        tic;    
        
        % Set log level (integer parameter)
        %param.MSK_IPAR_LOG = 1;
        % Select interior-point optimizer... (integer parameter)
        param.MSK_DPAR_INTPNT_CO_TOL_REL_GAP = 1.0e-6;
        
        if QUIET
            [r,res] = mosekopt('minimize echo(0)',prob, param);
        else
            [r,res] = mosekopt('minimize',prob, param);
        end
        
        time_solve = toc;
        if  strcmp(res.sol.itr.prosta, 'PRIMAL_AND_DUAL_FEASIBLE')        
            cost = res.sol.itr.pobjval;
        else
            cost = NaN;
        end

        package = res;
        [x,y] = convert_mosek2sedumi_var(res.sol.itr,prob.bardim);
        
    else
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
    
    %check for SDP optimality
    
    if isnan(cost)
        sdp_opt = false;
    else
%         x_rec = decomposed_recover(x, info);
%         [sdp_opt, cone_valid] = check_sdp_opt(x_rec, y, model.A, model.b, model.c, model.K, cone, dual);
        sdp_opt = true
    end
    
    Hout = sqrt(cost);
end

