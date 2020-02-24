load('sos_quartic_200.mat', 'model');
modelD = struct;
model.c =  model.C;
model = rmfield(model, 'c');

 [modelD.A, modelD.b, modelD.C, modelD.K, info] = ...
     dd_star_convert(model.A,model.b,model.C, model.K);

% As = model.A(:, 3:end);
% N = model.K.s;
% Nd = N *(N+1)/2;
% 
% %Constraints
% A_proc =  triu_matrix(As, 1);
% 
% 
% 
% %Additional Constraints
% rays = dd_extreme_rays(Ksi);
% raysT  = triu_matrix(rays', 1);
% 
%         A_rel_free{PSDind} = rays_svec';
%         A_rel_lin{PSDind} = -speye(num_dd);
% 
% %Cone
% Knew.f = model.K.f + Ksi*(Ksi+1)/2;
% Knew.l = model.K.l + Ksi^2;
% Knew.q = [];
% Knew.s = [];
% 
% %Cost
% C_proc = triu_vector(model.C, 1);
% 
% 
% modelD = struct;
% modelD.A = A_proc;
% modelD.C = C_proc;



%prob =  convert_sedumi2mosek(model.A,model.b,model.C, model.K);
%     dd_star_convert(model.A,model.b,model.c, model.K);
% [modelD.A, modelD.b, modelD.c, modelD.K, info] = ...
%     dd_star_convert(model.A,model.b,model.c, model.K);
%  
%[x, y, info_psd] = sedumi(model.A, model.b, model.C, model.K);
%s_rec = model.C - model.A'*y;


%  
 pars.fid = 0;
[xd, yd, info_dds] = sedumi(modelD.A, modelD.b, modelD.C, modelD.K, pars);

[Mi, L] = tri_indexer(model.K.s(1));
xs = xd(3:16);
x_ret = xs(Mi);
reshape(x_ret, 4, 4)
