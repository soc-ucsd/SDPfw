%load('sea_star_H2_small.mat')
%load('sea_star_Hinf0_small.mat')
%load('sea_star_Hinf0_medium.mat')

%fname = 'sea_star_H2_tiny.mat';
%fname = 'sea_star_Hinf0_tiny.mat';
%fname = 'sea_star_H2_small.mat';
%fname = 'sea_star_Hinf0_medium.mat';
fname = 'sea_star_Hinf0_wide_small.mat';
%fname = 'sea_star_Hinf0_large.mat';
%fname = 'sea_star_Hinf0_verylarge.mat';
%fname = 'sea_star_Hinf0_giant.mat';
%fname = 'sea_star_Hinf0_wide_small.mat';
[filepath,name,ext] = fileparts(fname);

outname = strcat(filepath,'output_',name,'_LOWER',ext);
load(fname, 'model_upper');
model = model_upper;


%running options
RUN_PSD = 1;
PRIMAL = 1;
DUAL = 1;

%medium
%thresh = [0, 11, 35, 100];
%cones = {'dd', 1, 3, 5, 9, 18, 36};

thresh = [0, 11, 35, 100];
cones = {'dd', 5, 10};

%thresh = [0,11];
%cones = {'dd', 2, 4, 6};
%thresh = [0];
%cones =  {'dd'};
%large

%thresh = [0, 11, 60, 100];
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
%cones = {'dd', 2, 4, 6};
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
QUIET = 1;

output = NaN*ones(Ncones, Nthresh);

Ks = model.K.s;

CONE = cell(Ncones, Nthresh);
RES  = cell(Ncones, Nthresh);
CONE_dual = cell(Ncones, Nthresh);
%RES_dual = cell(Ncones, Nthresh);

for i = 1:Ncones
    for j = 1:Nthresh
        CONE{i,j} = struct;
        CONE{i,j}.cone = cone_list(Ks, thresh(j), cones{i});
    end
end

CONE_dual = CONE; 

cone = cell(length(Ks), 1);

save(outname, 'cones', 'thresh');

for i = 1:Ncones
    for j = 1:Nthresh
        if PRIMAL
            %upper bound
            [CONE{i,j}.Hout, RES{i,j}, CONE{i,j}.time_solve, CONE{i,j}.time_convert]...
                = run_model_star(model, CONE{i,j}.cone, 0, use_mosek, QUIET);              

            fprintf('Cone: %s \t Thresh:  %d \t Hinf upper: %0.3f \t \t Time Solve: %0.1f \t Time Convert: %0.1f\n', ...
                num2str(cones{i}), thresh (j), CONE{i,j}.Hout, CONE{i,j}.time_solve, CONE{i,j}.time_convert)
        end
        
        if DUAL
            %lower bound?
            [CONE_dual{i,j}.Hout, RES{i,j}, CONE_dual{i,j}.time_solve, CONE_dual{i,j}.time_convert]...
                = run_model_star(model, CONE_dual{i,j}.cone, 1, use_mosek, QUIET);

            fprintf('Cone: %s \t Thresh:  %d \t Hinf lower: %0.3f \t \t Time Solve: %0.1f \t Time Convert: %0.1f\n', ...
                num2str(cones{i}), thresh (j), CONE_dual{i,j}.Hout, CONE_dual{i,j}.time_solve, CONE_dual{i,j}.time_convert)
            end
            save(outname,  '-append',  'CONE', 'CONE_dual');
        %save(outname, 'CONE', 'CONE0', 'cones', 'thresh')

    end
end



if PRIMAL
    cost = cellfun(@(x) x.Hout, CONE);
    time = cellfun(@(x) x.time_solve + x.time_convert, CONE);
    save(outname, '-append', 'cost', 'time')
end

if DUAL
    cost_dual = cellfun(@(x) x.Hout, CONE_dual);
    time_dual = cellfun(@(x) x.time_solve + x.time_convert, CONE_dual);
    save(outname, '-append', 'cost_dual', 'time_dual')
end
if RUN_PSD
    [CONE0.Hout, RES0, CONE0.time_solve, CONE0.time_convert]...
            = run_model_star(model, 'psd', use_mosek);
    fprintf('Cone: PSD \t Thresh: N\\A \t Hinf: %3f\t \t Time Solve: %0.1f \t Time Convert: %0.1f\n', ...
        CONE0.Hout,CONE0.time_solve, CONE0.time_convert )
    save(outname, '-append', 'CONE0')
end