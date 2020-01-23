%LR test 
sos= 1;
if sos
    load('LR_120_sos.mat')
    outname = 'LR120_output_uncons_sos.mat';
    thresh = [0, 100];
    cones = {'dd', 'sdd', 2, 3, 5, 11, 21, 40, 'psd'};
    model = model_csp;
    model.c = -model.C;
    support_LR(model, outname, cones, thresh);

else
    load('LR_120.mat')
    ones = {'dd'};
    cones = {'dd', 'sdd', 2, 3, 5, 6, 11, 20, 30, 40, 'psd'};

    thresh = [0, 11, 45, 100];
    outname_unc = 'LR120_output_uncons.mat';
    support_LR(model_unc, outname_unc, cones, thresh);


    outname_c = 'LR120_output_cons.mat';
    support_LR(model_c, outname_c, cones, thresh);

end



%cones = {20};

%thresh = 100;



function support_LR(model, outname, cones, thresh)
    Ncones = length(cones);
    Nthresh = length(thresh);

    CONE = cell(Ncones, Nthresh);
    RES  = cell(Ncones, Nthresh);

    for i = 1:Ncones
        for j = 1:Nthresh
            CONE{i,j} = struct;
            CONE{i,j}.cone = cone_list(model.K.s, thresh(j), cones{i});
        end
    end

    use_mosek = 1;

    cost = NaN*ones(Ncones, Nthresh);

    for i = 1:Ncones
        for j = 1:Nthresh
            [CONE{i,j}.cost, RES{i,j}, CONE{i,j}.time_solve, CONE{i,j}.time_convert]...
                = run_model_LR(model, CONE{i,j}.cone, use_mosek);              

            cost(i, j) = CONE{i,j}.cost;
                  %output
                 %{cones{i}, thresh(j), output(i,j)}
             fprintf('Cone: %s \t Thresh:  %d \t Cost: %0.3f \t \t Time Solve: %0.1f \t Time Convert: %0.1f\n', ...
                 num2str(cones{i}), thresh (j), cost(i,j), CONE{i,j}.time_solve, CONE{i,j}.time_convert)
    %         else
    %             CONE{i,j}.Hout = NaN;        
    %         end
    save(outname, 'CONE', 'cones', 'thresh')
        end
    end
end
% 
% CONE0 = struct;
% RES0 = struct;
% 
% [CONE0.Hout, RES0, CONE0.time_solve, CONE0.time_convert]...
%             = run_model_star(model, 'psd', use_mosek);
% fprintf('Cone: PSD \t Hinf: %3f\n', CONE0.Hout)
% save(outname, 'CONE', 'CONE0', 'cones', 'thresh')
