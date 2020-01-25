%fname = 'LR12_trans_box.json';
fname = 'LR120_trans_box.json';
%fname = 'LR120_cons_1.json';

%This is the golden file.
%Protect it at all costs
%fname = 'LR120_cons_3.json';


data = jsondecode(fileread(fname));

N = data.n;

%N_f = length(data.fbasis)/N;
fbasis = sparse(reshape(data.fbasis, N, []));
gbasis = cell(N, 1);
for i = 1:N
    gb_curr = data.gbasis(i, :);
    gbasis{i} = sparse(reshape(gb_curr, N, []));
end


fblocks = data.fblocks;
gblocks = data.gblocks;
% 
% %formulas 22 and 23 in the TSSOS paper in order to get the Psatz in YALMIP
% %each block in fblocks indexes into the support set in fbasis
% %likewise for gblocks
% 
% 
 x = sdpvar(N, 1);
 lower = sdpvar(1);
%gb = mon_basis(x, gbasis{1})
%fb = mon_basis(x, fbasis);

%start with the Psatz expression
%f_ref = f - lower;

f_psatz = 0;
%F = [];
%s = [];
%c = [];
Q = [];
F = [];
% 
suppf = reshape(data.ssupp{1}, N, []);
f_basis = mon_basis(x, suppf);
coef  = data.coe{1};
f = coef'*f_basis;

f_ref = f - lower;

%unconstrained terms
for i = 1:length(fblocks)
    x_basis = mon_basis(x, fbasis(:, fblocks{i}));
    Qi = sdpvar(length(x_basis));
    F = [F; Qi >= 0];    
    pi = x_basis'*Qi*x_basis;
    f_psatz = f_psatz + pi;
    %[si,ci] = polynomial(x_basis,2);
    %s = [s; si];
    %c = [c; ci];
    %F = [F; sos(si)];
end

%constrained terms
for j = 1:length(gblocks)
    gbasis_curr = gbasis{j};
    gblocks_curr = gblocks{j};
    pgj = 0;
    
    suppj = reshape(data.ssupp{j+1}, N, []);
    coej  = data.coe{j+1};
    g_basis = mon_basis(x, suppj);
    gj = coej'*g_basis;
    
    for i = 1:length(gblocks_curr)
        x_basis = mon_basis(x, gbasis_curr(:, gblocks_curr{i}));
        Qi = sdpvar(length(x_basis));
        F = [F; Qi >= 0];    
        pi = x_basis'*Qi*x_basis;
        %accumulated constrained multiplier
        pgj = pgj + pi;
        %f_psatz = f_psatz + pi;
    end
    
    f_psatz = f_psatz + pgj * gj;
end

F = [F; coefficients(f_ref - f_psatz, x) == 0];

obj = -lower;

opts = sdpsettings('solver', 'sedumi');

%sol = optimize(
[model, recoverymodel] = export(F, obj, opts);
%optimize(F, obj, opts);
[~,name, ~] = fileparts(fname);
outname = strcat(name, "_sos.mat");
save(outname, "model", "data")