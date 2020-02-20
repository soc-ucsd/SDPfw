%load('sea_star_H2_small.mat')
%load('sea_star_Hinf0_small.mat')
%load('sea_star_Hinf0_medium.mat')

LOWER = 1;
CONVERT = 0;
%fname = 'sea_star_H2_tiny.mat';
%fname = 'sea_star_Hinf0_tiny.mat';
%fname = 'sea_star_H2_small.mat';
%fname = 'sea_star_Hinf0_medium.mat';
%fname = 'sea_star_Hinf0_large.mat';
%fname = 'sea_star_Hinf0_verylarge.mat';
%fname = 'sea_star_Hinf0_giant.mat';
fname = 'sea_star_Hinf0_wide_small.mat';
[filepath,name,ext] = fileparts(fname);



if LOWER
    outname = strcat(filepath,'output_',name,'_LOWER',ext);
    load(fname, 'model_lower');
    model = model_lower;
else
    outname = strcat(filepath,'output_',name,ext);
    if CONVERT
%         load(fname, 'model_lower');
%         [F, h] = sedumi2yalmip(model_lower.A, model_lower.b,model_lower.c,model_lower.K);
%         [Fd, hd] = dualize(F, h, 0);
%         model_upper = export(Fd, -hd, sdpsettings('solver', 'sedumi'));
%         model_upper.c = model_upper.C;
%         model_upper_c = model_upper;
%         save(fname, '-append', 'model_upper_c');
%         outname = strcat(filepath,'output_CONVERT_',name,ext);
    else
        load(fname, 'model_upper');
    end
    
    model = model_upper;
end






% %medium
% thresh = [0, 11, 35, 100];
% cones = {'dd', 1, 3, 5, 9, 18, 36};

%large
%
thresh = [0, 11, 60, 100];
%thresh = [100];
%thresh = [0, 5, 11];
%thresh = [0, 11, 40];
%thresh = [0, 11, 41, 63];

%cones = {'dd', 1, 5, 7, 15, 30, 50};
%cones = {'dd', 1, 3, 5, 8, 15, 30, 55, 70};
%cones = {1, 3, 5, 8, 15, 30, 55, 70};
%cones = {73};
%cones = {'dd'};
%cones = {};
%cones = {'dd', 'sdd', 3, 6, 12, 18, Inf};
cones = {'dd', 2, 4, 6};
%cones = {1,2,3,4,5, 6};
%cones = {1,2, 4, 5, 8, 10};
%thresh = [0, 22, 50, Inf];
%thresh = [0, 11, 35, Inf];
%thresh = [0, 11, 30, 50,  Inf];
%thresh = [0, 11, 33,  56,  Inf];
%thresh = [0, 11, 33,  56];

Ncones = length(cones);
Nthresh = length(thresh);



use_mosek = 1;

output = NaN*ones(Ncones, Nthresh);
% % 
% if DUAL    
%     model_orig = struct;
%     model_orig.A = -LOP.A';
%     model_orig.b = -LOP.c;
%     model_orig.c = -LOP.b;
%     model_orig.K = LOP.J;
%     
%     [F, h] = sedumi2yalmip(model_orig.A, model_orig.b,model_orig.c,model_orig.K);
%     %[Fd, hd] = dualize(F, h)
%     [modelp,~] = export(F, h, sdpsettings('solver', 'sedumi', 'removeequalities', 2));
%     [Fp, hp] = sedumi2yalmip(modelp.A, modelp.b, modelp.C, modelp.K);
%     [Fd,hd,X,t,err] = dualize(F,h,0);
%     [model,~] = export(Fd, -hd, sdpsettings('solver', 'sedumi', 'removeequalities', 1));


%       LOP_dual = model;
%       save(fname, '-append', 'model_dual');
% end

%model.c = model.C;
% else
%     model = struct;
%     model.A = -LOP.A';
%     model.b = -LOP.c;
%     model.c = -LOP.b;
%     model.K = LOP.J;
% end


%Js = model.J.s';
Js = model.K.s;

CONE = cell(Ncones, Nthresh);
RES  = cell(Ncones, Nthresh);
CONE0 = struct;
RES0 = struct;
for i = 1:Ncones
    for j = 1:Nthresh
        CONE{i,j} = struct;
        CONE{i,j}.cone = cone_list(Js, thresh(j), cones{i});
    end
end

cone = cell(length(Js), 1);

%[A, b, c, K, ~] = decomposed_subset(-model.A',-model.c,-model.b,model.J, cone);

for i = 1:Ncones
    for j = 1:Nthresh
        [CONE{i,j}.Hout, RES{i,j}, CONE{i,j}.time_solve, CONE{i,j}.time_convert]...
            = run_model_star(model, CONE{i,j}.cone, use_mosek);              
             
        
        if LOWER
            CONE{i,j}.Hout = imag(CONE{i,j}.Hout);
        end
        
        output(i, j) = CONE{i,j}.Hout;
              %output
             %{cones{i}, thresh(j), output(i,j)}
         fprintf('Cone: %s \t Thresh:  %d \t Hinf: %0.3f \t \t Time Solve: %0.1f \t Time Convert: %0.1f\n', ...
             num2str(cones{i}), thresh (j), output(i,j), CONE{i,j}.time_solve, CONE{i,j}.time_convert)
%         else
%             CONE{i,j}.Hout = NaN;        
%         end
save(outname, 'CONE', 'CONE0', 'cones', 'thresh')
    end
end
% 



cost = cellfun(@(x) x.Hout, CONE);
time = cellfun(@(x) x.time_solve + x.time_convert, CONE);


[CONE0.Hout, RES0, CONE0.time_solve, CONE0.time_convert]...
            = run_model_star(model, 'psd', use_mosek);

if LOWER
    CONE0.Hout = imag(CONE0.Hout);
end
        fprintf('Cone: PSD \t Hinf: %3f\n', CONE0.Hout)
save(outname, 'CONE', 'CONE0', 'cones', 'thresh')





% CONE73=struct;
% cone73 = cone_list(Js, 100, 73);
% [CONE73.Hout, RES73, CONE73.time_solve, CONE73.time_convert]...
%             = run_model_star(model, cone73, use_mosek);


