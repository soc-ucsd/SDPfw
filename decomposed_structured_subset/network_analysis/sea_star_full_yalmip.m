%load('sea_star_H2_tiny.mat', 'n', 'm','d','Sys')
%load('sea_star_Hinf0_small.mat', 'n', 'm','d','Sys')
fname = 'sea_star_Hinf0_verylarge.mat';
%fname = 'sea_star_H2_tiny.mat';
[filepath,name,ext] = fileparts(fname);
outname = strcat(filepath,'output_dense_',name,ext);

load(fname, 'n', 'm','d','Sys');

cones = {'dd', 1, 3, 5, 8, 15, 30, 55, 70};
%cones = {'dd', 10};
CONE = cell(length(cones), 1);
for k = 1:length(cones)
    cone = cones(k);
    CONE{k}.cone = cone;
    
    tic
    epsilon = 0.01;
    P = [];
    for i = 1:length(n)
        P = blkdiag(P,sdpvar(n(i)));
    end
    P_reg = P - epsilon*eye(sum(n));
    N = sum(n);

    num_input  = size(Sys.globalB, 2);
    num_output = size(Sys.globalC, 1);
    num_state  = size(Sys.globalA, 2);

    Sys.globalD = sparse(num_output, num_input);

    gamma2 = sdpvar(1);
    Bounded_Real = -[P*Sys.globalA+Sys.globalA'*P + Sys.globalC'*Sys.globalC, P*Sys.globalB + Sys.globalC'*Sys.globalD; 
                        Sys.globalB'*P + Sys.globalD'*Sys.globalC, Sys.globalD'*Sys.globalD-gamma2*eye(sum(m))]  - epsilon*eye(sum(n)+sum(m));
    NB = length(Bounded_Real);

    %cone = 'dd';                
    %cone = 10;                
    Constraint1 = [];
    Constraint2 = [];

    if strcmp(cone, 'psd')
        Constraint1 = [P_reg >= 0];
        Constraint2 = [ Bounded_Real >= 0];
    elseif strcmp(cone, 'dd')
        for i = 1:N
            i_curr = (1:N);
            i_curr(i) = [];
            Constraint1 = [Constraint1, P_reg(i, i) >= sum(abs(P_reg(i_curr, i)))];
        end

        for i = 1:NB
            i_curr = (1:NB);
            i_curr(i) = [];
            Constraint2 = [Constraint2, Bounded_Real(i, i) >= sum(abs(Bounded_Real(i_curr, i)))];
        end
    elseif isnumeric(cone)  && (cone > 0)      

        [Var1, Constraint1] = factorwidth_yalmip(P_reg, cone, 0);
        [Var2, Constraint2] = factorwidth_yalmip(Bounded_Real, cone, 0);

    end

    Constraint = [Constraint1, Constraint2];
    Cost = gamma2;
    CONE{k}.time_convert = toc;
    tic
    opts = sdpsettings('solver', 'mosek', 'verbose', 1);
    %opts = sdpsettings('solver', 'mosek', 'verbose', 0);
    sol = optimize(Constraint, Cost, opts);

    gamma = sqrt(value(gamma2));
    CONE{k}.time_solve = toc;
    CONE{k}.Hout = gamma;
    fprintf('Cone: %s \t Hinf: %0.3f \t \t Time Solve: %0.1f \t Time Convert: %0.1f\n', ...
        num2str(cone), CONE{k}.Hout, CONE{k}.time_solve, CONE{k}.time_convert)
    save(outname, 'CONE', 'cones');
end

function [Var, Constraint] = factorwidth_yalmip(M, cone, dual)
    %factor width and duals for yalmip
    nop = max(floor(length(M)/cone), 1);
    SizeU = ceil(length(M)/nop);
    SizeL = floor(length(M)/nop);
    if SizeU == SizeL
        alpha = ones(nop,1)*SizeU;
    else
        x = (length(M) - SizeL*nop)./(SizeU-SizeL);
        alpha = [ones(x,1)*SizeU;ones(nop-x,1)*SizeL];
    end
    alpha_sum = cumsum(alpha);
    alpha_sum_low= [0; alpha_sum(1:end-1)]+1;
    
    if length(alpha) > 1
        Num = nchoosek(length(alpha),2);    

        Constraint=[];
        if dual            
            Var = {};
        else        
            %Accum = zeros(length(M));
            %Constraint = [Accum == M];   
            Var = cell(Num+1, 1);        
        end

        k = 1;
        for i = 1:length(alpha)
            for j = (i+1):length(alpha)          
                block1 = alpha_sum_low(i):alpha_sum(i);
                block2 = alpha_sum_low(j):alpha_sum(j);
                block_joint = [block1, block2];
                numvar = length(block_joint);
                
                
                
                if dual
                    Constraint = [Constraint, M(block_joint, block_joint) >= 0];
                else
                    Var{k} = sdpvar(numvar, numvar);
                    E = sparse(1:length(block_joint), block_joint, ones(length(block_joint),1), length(block_joint),length(M));
                
                    if (i==1) && (j==2)
                        Accum = E'*Var{k}*E;
                    else
                        Accum = Accum + E'*Var{k}*E;
                    end
                    %Accum(block_joint, block_joint) = Accum(block_joint, block_joint)  + Var{k};                
                    Constraint = [Constraint, Var{k} >=0];
                end



                k = k + 1;

            end
        end

        if ~dual
            Constraint= [Constraint, M==Accum];
        end
    else
        Constraint = [M >= 0];
        Var = {};
    end
end   