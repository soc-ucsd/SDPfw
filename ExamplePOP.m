
%% Example 1: Lower bounds on polynomial optimization problems 
% Broyden tridiagonal function 


Num   = [10,15,20,30,40,50];                % number of variables
Lower = zeros(length(Num),4);  % SOS, CSOS, SDSOS, DSOS
Time  = zeros(length(Num),4);

for ind = 6:length(Num)

%% Method 1: solution via Yalmip - orginal SOS
    N = Num(ind);
    
    fprintf('Number of dimension: %d \n',N)
    
    fprintf('Testing original SOS ... : ')
    d = 2;
    x = sdpvar(N,1);
    f = ((3-2*x(1))*x(1)-2*x(2)+1)^2;
    for i =2:N-1
        f = f + ((3-2*x(i))*x(i)-x(i-1)-2*x(i+1)+1)^2;
    end
    f = f+ ((3-2*x(N))*x(N)-x(N-1)+1)^2;
    
    gamma = sdpvar(1);
    
    F = sos(f-gamma*(x'*x));
    %opts = sdpsettings('solver','mosek');
    opts = sdpsettings('solver','sedumi');
    opts.sos.newton = 0;
    opts.sos.congruence = 0;
    opts.sos.csp = 0; 
    
    opts.Num = N;
    
    %if N < 30
    [sol,v,Q] = solvesos(F,-gamma,opts);
    %end
    
    fprintf('finished  %4.2f \n', sol.solvertime);
    
%% Method 2: solution via Yalmip - decomposing SOS
    fprintf('Testing decomposing SOS ... : ')
    gamma1 = sdpvar(1);
    F = sos(f-gamma1*(x'*x));
    opts.sos.csp = 1;                    % correlative sparsity pattern
    [sol1,v,Q] = solvesos(F,-gamma1,opts);
    fprintf('finished  %4.2f \n', sol1.solvertime);
    
%% Method 3: solution via Spotless - DSOS/SDSOS
    y      = msspoly('y',N);
    g = ((3-2*y(1))*y(1)-2*y(2)+1)^2;
    for i =2:N-1
        g = g + ((3-2*y(i))*y(i)-y(i-1)-2*y(i+1)+1)^2;
    end
    g = g+ ((3-2*y(N))*y(N)-y(N-1)+1)^2;
    
    prog = spotsosprog;
    prog = prog.withIndeterminate(y);

    % New free variable gamma
    [prog,gamma2] = prog.newFree(1);
    % DSOS constraint
    prog = prog.withDSOS(g - gamma2*(y'*y)); % Only line that changes between DSOS,SDSOS,SOS programs

    % MOSEK options
    options = spot_sdp_default_options();
    options.solveroptions.verbose = 1;
    options.verbose = 1;
    options.solveroptions.MSK_IPAR_BI_CLEAN_OPTIMIZER = 'MSK_OPTIMIZER_INTPNT'; % Use just the interior point algorithm to clean up
    options.solveroptions.MSK_IPAR_INTPNT_BASIS = 'MSK_BI_NEVER'; % Don't use basis identification (it's slow)

    % Solve program
    sol2 = prog.minimize(-gamma2, @spot_mosek, options);
    if strcmp(sol2.status,'STATUS_PRIMAL_INFEASIBLE')
        opt_dsos = NaN;
    else
        opt_dsos = double(sol2.eval(gamma2));
    end
    %% SDSOS
    prog = spotsosprog;
    prog = prog.withIndeterminate(y);
    [prog,gamma3] = prog.newFree(1);
    % DSOS constraint
    prog = prog.withSDSOS(g - gamma3*(y'*y)); % Only line that changes between DSOS,SDSOS,SOS programs
    % Solve program
    sol3 = prog.minimize(-gamma3, @spot_mosek, options);
    opt_sdsos = double(sol3.eval(gamma3));
    
    %% statistics 
    Lower(ind,:) = [value(gamma),value(gamma1),opt_dsos,opt_sdsos];
    Time(ind,:)  = [sol.solvertime,sol1.solvertime,sol2.info.wtime,sol3.info.wtime];
    
    save POPbt

end