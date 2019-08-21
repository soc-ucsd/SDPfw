
%% Example 1: Lower bounds on polynomial optimization problems 
% Broyden tridiagonal function 


Num   = [10,15,20,25,30,35,40];                % number of variables
Lower = zeros(length(Num),4);  % SOS, CSOS, SDSOS, DSOS
Time  = zeros(length(Num),4);

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
    %opts = sdpsettings('solver','mosek');
    opts = sdpsettings('solver','sedumi');
    opts.sos.newton = 0;
    opts.sos.congruence = 0;
    opts.sos.csp = 0; 
    
    %opts.Num = N;
    
    [model,recoverymodel,diagnostic,internalmodel] = export(F,-gamma,opts);
    
    A = model.A;
    b = model.b;
    c = model.C;
    K = model.K;
    save([folder,'SedumiDataEx',num2str(N)],'A','b','c','K') 
    %if N < 30
    %[sol,v,Q] = solvesos(F,-gamma,opts);
    %end
    
    fprintf('finished  %4.2f \n', N);
    
end