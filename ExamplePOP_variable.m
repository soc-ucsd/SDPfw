%% Example 1: Lower bounds on polynomial optimization problems 
% Broyden tridiagonal function 

var_par = 3:10; % partition in the variables x
Num     = [10,15,20,25,30];                % number of variables ,35,40
Cost   = zeros(length(Num),length(var_par));  % SOS, CSOS, SDSOS, DSOS
Time    = zeros(length(Num),length(var_par));
Tcon  = zeros(length(Num),length(var_par));

folder = 'SeDuMiData\';

for ind = 1:length(Num)

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
    f = f + sum(x).^2;
    
    gamma = sdpvar(1);
    
    F = sos(f-gamma);
    
    fprintf('finished  %4.2f \n', N);
%% Method 2: solution via Yalmip
    for jj = 1:length(var_par)
        m = var_par(jj);
        if N > 2*m
            cliques = num2cell((0:m:(N-m))' + (1:m),2); %
            if (cliques{end}(end) < N)
                cliques{length(cliques)+1} = (cliques{end}(1)+m):N; 
            end

            if length(cliques) > 2
                combos  = combntns(1:length(cliques),2);
                I      = num2cell(combos,2);
                degree = 2;
                opts = sdpsettings('solver','mosek');
                opts.sos.newton = 0;
                opts.sos.congruence = 0;
                opts.sos.csp = 0;      
                % I am not sure if this monomial decomposition is correct.
                % Can you double check? 
                monoms = cellfun(@(C)[monolist(x(cliques{C(1)}),(degree+2)/2);monolist(x(cliques{C(2)}),(degree+2)/2)], I, 'UniformOutput', 0);
                try
                    [sol1,v1,Q1]   = solvesos(F,-gamma,opts,[gamma],{monoms});
                    Time(ind,jj) = sol1.solvertime;
                    Tcon(ind,jj) = sol1.yalmiptime;
                    Cost(ind,jj) = -value(gamma);
                catch
                    Time(ind,jj) = NaN;
                    Tcon(ind,jj) = NaN;
                    Cost(ind,jj) = NaN;
                end
            end
        end
        save([folder, 'Result_var_decomp'],'Time','Tcon','Cost')
    end

    yalmip('clear')

end