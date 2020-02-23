 function Sys = GenerateDynamics(Gp,Gc,n,m,d,Flag,Q,R)
% Generate a network of linear hetergoeneous systems
% Gp: Plant graph
% Gc: Communication graph
% n:  state dimension
% m:  input dimension
% d:  output dimension
% Flag: 0 deterministic dynamics
% Flag: 1 stability; 2 H2; 3 Hinf; 
%       4 stabilization, 5 H2 control; 6 Hinf control

if nargin <= 4
    d = [];
    Flag = 1;
end

[N,~] = size(Gp);
Sys.A       = cell(N,N);
Sys.globalA = zeros(sum(n),sum(n)); 
cumN        = cumsum([1;n(:)]);
cumM        = cumsum([1;m(:)]);
cumD        = cumsum([1;d(:)]);

if Flag == 0
    % Dynamic matrices
    n     = 2*ones(1,N);
    m     = ones(1,N);
    d     = n+m;
    cumN  = cumsum([1;n(:)]);
    cumM  = cumsum([1;m(:)]);
    cumD  = cumsum([1;d(:)]);
    Sys.A = cell(N,N);   % matrices for A
    Sys.B1 = cell(N,1);   % matrices for B
    Sys.B2 = cell(N,1);   % matrices for B
    Sys.C = cell(N,1);   % matrices for B
    Sys.D = cell(N,1);   % matrices for B
    
    Sys.globalA = zeros(sum(n),sum(n)); 
    Sys.globalB1 = zeros(sum(n),sum(m)); 
    Sys.globalB2 = zeros(sum(n),sum(m)); 
    Sys.globalC = zeros(sum(d),sum(n)); 
    Sys.globalD = zeros(sum(d),sum(m));     

    % Dynamical part
    for i = 1:N
        Sys.B1{i} = [0;1];
        Sys.B2{i} = [0;1];
        Sys.C{i} = [Q^(1/2); zeros(m(i),n(i))];
        Sys.D{i} = [zeros(n(i),m(i)); R^(1/2)];
        Sys.A{i,i} = [1 1; 1 2];
        Sys.globalA(cumN(i):cumN(i+1)-1,cumN(i):cumN(i+1)-1) = Sys.A{i,i};
        Sys.globalB1(cumN(i):cumN(i+1)-1,cumM(i):cumM(i+1)-1) = Sys.B1{i};
        Sys.globalB2(cumN(i):cumN(i+1)-1,cumM(i):cumM(i+1)-1) = Sys.B2{i};
        Sys.globalC(cumD(i):cumD(i+1)-1,cumN(i):cumN(i+1)-1) = Sys.C{i};
        Sys.globalD(cumD(i):cumD(i+1)-1,cumM(i):cumM(i+1)-1) = Sys.D{i};
        for j = 1:N
            if Gp(i,j) == 1 && i~=j
                Sys.A{i,j} = exp(-(i-j)^2/10)*eye(n(i));     %% node j has influence on node i
                Sys.globalA(cumN(i):cumN(i+1)-1,cumN(j):cumN(j+1)-1) = Sys.A{i,j};
            end   
        end
    end
end

if Flag == 1          % stabilization test
    for i = 1:N
       Sys.A{i,i} = rand(n(i));
       Sys.globalA(cumN(i):cumN(i+1)-1,cumN(i):cumN(i+1)-1) = Sys.A{i,i};
        for j = 1:N
           if Gp(i,j) == 1 && i ~= j
               Sys.A{i,j} = rand(n(i),n(j));
               Sys.globalA(cumN(i):cumN(i+1)-1,cumN(j):cumN(j+1)-1) = Sys.A{i,j};
           end
        end
    end
    % guanrantee stability
    maxLambda = max(eig(Sys.globalA));
    Sys.globalA = Sys.globalA - (maxLambda + 5)*eye(sum(n));
    for i = 1:N
        Sys.A{i,i} = Sys.A{i,i} - (maxLambda + 5 )*eye(n(i));
    end 
elseif Flag == 2  % H2 test
    Sys.B       = cell(N,1);
    Sys.C       = cell(N,1);
    Sys.globalB = zeros(sum(n),sum(m)); 
    Sys.globalC = zeros(sum(d),sum(n));
    for i = 1:N
       Sys.A{i,i} = rand(n(i));
       Sys.B{i}   = rand(n(i),m(i));
       Sys.C{i}   = rand(d(i),n(i));
       Sys.globalA(cumN(i):cumN(i+1)-1,cumN(i):cumN(i+1)-1) = Sys.A{i,i};
       Sys.globalB(cumN(i):cumN(i+1)-1,cumM(i):cumM(i+1)-1) = Sys.B{i};
       Sys.globalC(cumD(i):cumD(i+1)-1,cumN(i):cumN(i+1)-1) = Sys.C{i};
        for j = 1:N
           if Gp(i,j) == 1 && i ~= j
               Sys.A{i,j} = rand(n(i),n(j));
               Sys.globalA(cumN(i):cumN(i+1)-1,cumN(j):cumN(j+1)-1) = Sys.A{i,j};
           end
        end
    end
    % guanrantee stability
    maxLambda = max(eig(Sys.globalA));
    Sys.globalA = Sys.globalA - (maxLambda + 1)*eye(sum(n));
    for i = 1:N
        Sys.A{i,i} = Sys.A{i,i} - (maxLambda + 1 )*eye(n(i));
    end 
elseif Flag == 3  % Hinf test
    Sys.B       = cell(N,1);
    Sys.C       = cell(N,1);
    Sys.D       = cell(N,1);
    Sys.globalB = zeros(sum(n),sum(m)); 
    Sys.globalC = zeros(sum(d),sum(n));
    Sys.globalD = zeros(sum(d),sum(m));
    for i = 1:N
       Sys.A{i,i} = rand(n(i));
       Sys.B{i}   = rand(n(i),m(i));
       Sys.C{i}   = rand(d(i),n(i));
       Sys.D{i}   = rand(d(i),m(i));
       Sys.globalA(cumN(i):cumN(i+1)-1,cumN(i):cumN(i+1)-1) = Sys.A{i,i};
       Sys.globalB(cumN(i):cumN(i+1)-1,cumM(i):cumM(i+1)-1) = Sys.B{i};
       Sys.globalC(cumD(i):cumD(i+1)-1,cumN(i):cumN(i+1)-1) = Sys.C{i};
       Sys.globalD(cumD(i):cumD(i+1)-1,cumM(i):cumM(i+1)-1) = Sys.D{i};
        for j = 1:N
           if Gp(i,j) == 1 && i ~= j
               Sys.A{i,j} = rand(n(i),n(j));
               Sys.globalA(cumN(i):cumN(i+1)-1,cumN(j):cumN(j+1)-1) = Sys.A{i,j};
           end
        end
    end
    % guanrantee stability
    maxLambda = max(eig(Sys.globalA));
    Sys.globalA = Sys.globalA - (maxLambda + 5)*eye(sum(n));
    for i = 1:N
        Sys.A{i,i} = Sys.A{i,i} - (maxLambda + 5 )*eye(n(i));
    end 
    
elseif Flag == 4  % stablization
    Sys.B       = cell(N,1);
    Sys.globalB = zeros(sum(n),sum(m));
    for i = 1:N
       Sys.A{i,i} = rand(n(i));
       Sys.B{i}   = rand(n(i),m(i));
       Sys.globalA(cumN(i):cumN(i+1)-1,cumN(i):cumN(i+1)-1)  = Sys.A{i,i};
       Sys.globalB(cumN(i):cumN(i+1)-1,cumM(i):cumM(i+1)-1) = Sys.B{i};
        for j = 1:N
           if Gp(i,j) == 1 && i ~= j
               Sys.A{i,j} = rand(n(i),n(j));
               Sys.globalA(cumN(i):cumN(i+1)-1,cumN(j):cumN(j+1)-1) = Sys.A{i,j};
           end
        end
    end   
    
elseif Flag == 5 % for H2 control  
    Sys.B1       = cell(N,1);
    Sys.B2       = cell(N,1);
    Sys.C        = cell(N,1);
    Sys.D        = cell(N,1);
    Sys.globalB1 = zeros(sum(n),sum(m)); 
    Sys.globalB2 = zeros(sum(n),sum(m));
    Sys.globalC  = zeros(sum(d),sum(n));
    Sys.globalD  = zeros(sum(d),sum(m));
    for i = 1:N
       Sys.A{i,i} = rand(n(i));
       Sys.B1{i}  = rand(n(i),m(i));
       Sys.B2{i}  = rand(n(i),m(i));
       Sys.C{i}   = rand(d(i),n(i));
       Sys.D{i}   = 2*rand(d(i),m(i));
       Sys.globalA(cumN(i):cumN(i+1)-1,cumN(i):cumN(i+1)-1)  = Sys.A{i,i};
       Sys.globalB1(cumN(i):cumN(i+1)-1,cumM(i):cumM(i+1)-1) = Sys.B1{i};
       Sys.globalB2(cumN(i):cumN(i+1)-1,cumM(i):cumM(i+1)-1) = Sys.B2{i};
       Sys.globalC(cumD(i):cumD(i+1)-1,cumN(i):cumN(i+1)-1)  = Sys.C{i};
       Sys.globalD(cumD(i):cumD(i+1)-1,cumM(i):cumM(i+1)-1)  = Sys.D{i};
        for j = 1:N
           if Gp(i,j) == 1 && i ~= j
               Sys.A{i,j} = rand(n(i),n(j));
               Sys.globalA(cumN(i):cumN(i+1)-1,cumN(j):cumN(j+1)-1) = Sys.A{i,j};
           end
        end
    end   
elseif Flag == 6
    
end






end

