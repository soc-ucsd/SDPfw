%load('sea_star_H2_small.mat')
%load('sea_star_Hinf0_small.mat')
%load('sea_star_Hinf0_medium.mat')

%fname = 'sea_star_H2_tiny.mat';
%fname = 'sea_star_Hinf0_tiny.mat';
%fname = 'sea_star_H2_small.mat';
%fname = 'sea_star_Hinf0_medium.mat';
fname = 'sea_star_Hinf0_large.mat';
[filepath,name,ext] = fileparts(fname);
outname = strcat(filepath,'output_',name,ext);

load(fname);

Js = LOP.J.s';




% %medium
% thresh = [0, 11, 35, 100];
% cones = {'dd', 1, 3, 5, 9, 18, 36};

%large
%thresh = [0, 11, 40, 75, 100];
thresh = [0, 11, 41, 63];

cones = {'dd', 1, 3, 5, 7, 15, 30, 50};
%cones = {'dd', 1, 3, 5, 8, 15, 30, 60, 85};
%cones = {'dd'};

%cones = {'dd', 'sdd', 3, 6, 12, 18, Inf};
%cones = {'dd', 6, 18, Inf};
%thresh = [0, 22, 50, Inf];
%thresh = [0, 11, 35, Inf];
%thresh = [0, 11, 30, 50,  Inf];
%thresh = [0, 11, 33,  56,  Inf];
%thresh = [0, 11, 33,  56];

Ncones = length(cones);
Nthresh = length(thresh);

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

use_mosek = 1;

output = NaN*ones(Ncones, Nthresh);





for i = 1:Ncones
    for j = 1:Nthresh
        [CONE{i,j}.Hout, RES{i,j}, CONE{i,j}.time_solve, CONE{i,j}.time_convert]...
            = run_model_star(LOP, CONE{i,j}.cone, use_mosek);              
             
        output(i, j) = CONE{i,j}.Hout;
              %output
             %{cones{i}, thresh(j), output(i,j)}
         fprintf('Cone: %s \t Thresh:  %d \t Hinf: %3f\n', num2str(cones{i}), thresh (j), output(i,j))
%         else
%             CONE{i,j}.Hout = NaN;        
%         end
save(outname, 'CONE', 'CONE0', 'cones', 'thresh')
    end
end
% 
fprintf('Cone: PSD \t Hinf: %3f\n', CONE0.Hout)
[CONE0.Hout, RES0, CONE0.time_solve, CONE0.time_convert]...
            = run_model_star(LOP, 'psd', use_mosek);
save(outname, 'CONE', 'CONE0', 'cones', 'thresh')
% [output ones(Ncones, 1)*CONE0.Hout]
% 






