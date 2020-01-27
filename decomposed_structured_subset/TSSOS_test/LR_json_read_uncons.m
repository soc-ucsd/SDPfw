%fname = 'LR12_trans_box.json';
%fname = 'LR120_trans_box.json';
fname = 'LR120_uncons_3.json';

%This is the golden file.
%Protect it at all costs
%fname = 'LR120_cons_3.json';
%fname = 'LR120_uncons_3.json';

data = jsondecode(fileread(fname));

disp("Read JSON")

N = data.n;

%N_f = length(data.basis)/N;
basis = sparse(reshape(data.basis, N, []));


blocks = data.blocks;

% 
% %formulas 22 and 23 in the TSSOS paper in order to get the Psatz in YALMIP
% %each block in blocks indexes into the support set in basis
% %likewise for gblocks
% 
% 
 x = sdpvar(N, 1);
 lower = sdpvar(1);
%gb = mon_basis(x, gbasis{1})
%fb = mon_basis(x, basis);

%start with the Psatz expression
%f_ref = f - lower;

f_psatz = 0;
%F = [];
%s = [];
%c = [];
Q = [];
F = [];
% 
disp("Starting with f")
supp = reshape(data.supp, N, []);
f_basis = mon_basis(x, supp);
coef  = data.coe;
f = coef'*f_basis;

f_ref = f - lower;

%unconstrained terms
for i = 1:length(blocks)
    if mod(i, 100)==0
        disp(strcat("starting f", num2str(i)))
    end
    x_basis = mon_basis(x, basis(:, blocks{i}));
    Qi = sdpvar(length(x_basis));
    F = [F; Qi >= 0];    
    pi = x_basis'*Qi*x_basis;
    f_psatz = f_psatz + pi;
    
    %[si,ci] = polynomial(x_basis,2);
    %s = [s; si];
    %c = [c; ci];
    %F = [F; sos(si)];
end

disp("finished with f")

disp("About to take Coefficients")
[cp] = coefficients(f_ref - f_psatz, x);
%F_lmi = F;
F = [F; cp == 0];

obj = -lower;

opts = sdpsettings('solver', 'mosek', 'savedebug', 1, 'savesolverinput', 1);
disp("Ready to Optimize")
diagnostics = optimize(F, obj, opts);
disp("Done Optimizing")
model = diagnostics.solverinput;
save("LR120_uncons_mosek_input.mat", "model")
%[model, recoverymodel] = export(F, obj, opts);
%optimize(F, obj, opts);
%[~,name, ~] = fileparts(fname);
%outname = strcat(name, "_sos.mat");
%save(outname, "model", "data")



