%% preliminary stuff

clear; clc; %close all
%define vec function
vec = @(x) x(:);
%add paths you need, excluding minConf - https://www.cs.ubc.ca/~schmidtm/Software/minConf.html
addpath(genpath('../pSPOT')); %spot operators package for faster linear-algebra in matlab
addpath(genpath('../spot')); %https://www.cs.ubc.ca/labs/scl/spot/

%% pick grid parameters
params.ns = 64; %sources
params.resolution = 5; %grid resolution -> set at 5km by now
%pick receiver grid
params.Gx = [70,165];
params.Gy = [70, 165]; 
%%%%%%%%%
params.Rx = (params.Gx(2)-params.Gx(1))/params.resolution+1;
params.Ry = (params.Gy(2)-params.Gx(1))/params.resolution+1; 
%up here, you need to store all the indices together so that things are
%ordered correctly
params.model_cutoff = 400+1; %build in buffer -> only used when manually binning data
%Note that if you want to interpolate real data, then you have to pick some
%perc not zero since then you can't compare recovered data
%Default for model data is 0.
params.perc = .2; %note: perc = 0 => see everything. perc = 1 => see nothing
params.figs_on = 1; %1 -> figs on, 0 -> figs off
%% Extract all the data, subsampling, and tensor sorting 
%do the subsampling part (not necessary for Ken's, unless you want to
%subsample observed data to get an  unseen_recov value)
params.interp = 'res'; %'tt' or 'res'
params.sub_used = 'ken'; %'ken' or 'rajiv' or 'real'
%subsamp matrix you want
params.norm_type = 'cb';%l2, l1, cb, lap note: lap still under construction
pre_processed_data = varproj_dataread(params); %extracts the data into a usable form
%change the data and get comparisons out (can comment out - used for
%checking efficacy of 2D projection operator)
% params.switchdata = 'n';  %yes to switch data with model projected data, no to keep regular data and just get out stats
% [pre_processed_data, opstats] = varproj_optest(params, pre_processed_data); 
%now do the subsampling step
sub_data = varproj_sub(params, pre_processed_data); 
 

%% Initialize and run the algorithm
%         params.sigma = 100; 
%         params.sigma =1e-10*norm(vec(b),1);
%         params.sigma = sqrt(.032^2*numel(sub_data.obs));
%         params.sigma = .5; 
% params.sigma = sqrt(0.0697^2*numel(sub_data.obs)); 
params.sigma = 0; 
params = varproj_setparms(params, sub_data, pre_processed_data); 

%% Run Algorithms
%(all take up the same number of iterations)
params.sup_out = 0; %1 -> do not display iteration progress when running the algorithm. 
interpolated_cb = varproj_alg(params, sub_data, vec); %with sigma = 0
%run FISTA (note options are for FISTA implementation)
options.iterations = params.stop_crit*params.iter_crit*10; 
options.verbosity = 2; %all
interpolated_fista = fista_nucnorm_lapreg(sub_data,params,options); % with sigma = 0

%run minFunc
params.iter_crit = 30; 
interpolated_lbfgs = minFnc_run(params, sub_data); 

params.iter_crit = 5;
params.stop_crit = 20; 
% params.sigma =  sqrt(0.0697^2*numel(sub_data.obs));
params.sigma = sqrt(.06^2*numel(sub_data.obs));
% params.sigma = 5.7843; 
% params.sigma = sqrt(.032^2*numel(sub_data.obs));
interpolated_cbnz = varproj_alg(params, sub_data, vec); 

%% comparison of only global and only local
% params.sigma = 5; 
params.norm_type = 'lap';
params = varproj_setparms(params, sub_data,pre_processed_data); 
interpolated_lap = varproj_alg(params, sub_data, vec); 

%low-rank alone is bad
params.norm_type = 'l2';
% params.sigma=1e-10*norm(vec(b),1);
params = varproj_setparms(params, sub_data,pre_processed_data); 
interpolated_l2 = varproj_alg(params, sub_data, vec); 

%% Do backend data visualization

params.figs_on = 1;
stats = varproj_datavis(params, sub_data, interpolated_cb, pre_processed_data, 'Var-Relax', 'vrs0'); 
stats_nz = varproj_datavis(params, sub_data, interpolated_cbnz, pre_processed_data, 'Var-Relax \sigma>0', 'vr');
stats_fist = varproj_datavis(params, sub_data, interpolated_fista, pre_processed_data, 'FISTA', 'fista'); 
stats_lbfgs = varproj_datavis(params, sub_data, interpolated_lbfgs, pre_processed_data, 'L-BFGS','lbfgs'); 
statslr = varproj_datavis(params, sub_data, interpolated_l2, pre_processed_data, 'LR Only', 'lr'); 
statslap = varproj_datavis(params, sub_data, interpolated_lap, pre_processed_data, 'Smooth Only', 'sm'); 


%% SVD and Objective value decay
varproj_figgen_svd(interpolated_cb.XM, 'VR - \sigma =0', ...
    interpolated_cbnz.XM, strcat('VR - \sigma = ', params.sigma), ...
    interpolated_fista.XM, 'FISTA', ...
    interpolated_lbfgs.XM, 'L-BFGS', ...
    interpolated_l2.XM, 'Low-Rank', ...
    interpolated_lap.XM, 'Smoothing', ...
    sub_data.real.res, 'Observed Res', ...
    sub_data.model.res, 'Model Res'); 
varproj_figgen_obj(interpolated_cb.obj_val, 'VR - \sigma =0', ...
    interpolated_cbnz.obj_val, strcat('VR - \sigma = ', params.sigma), ...
    interpolated_fista.obj_val, 'FISTA', ...
    interpolated_lbfgs.obj_val, 'L-BFGS', ...
    interpolated_l2.obj_val, 'Low-Rank', ...
    interpolated_lap.obj_val, 'Smoothing');

%% Bootstrapping

% params.bootruns = 1000; % however many bootstrapping runs you want
% params.figs_on = 0; 
% params.save_on = 1; 
% %run the bootstrapping function
% % [results, interpgrids, finalstats] = varproj_bootstrapping(params,pre_processed_data, vec); 
% %% FISTA and variable relaxation comparison
% %-doesn't work
% params.interp = 'res'; %'tt' or 'res'
% params.sub_used = 'ken'; %'ken' or 'rajiv' or 'real'
% %subsamp matrix you want
% params.norm_type = 'cb';%l2, l1, cb, lap note: lap still under construction
% params.sigma = 0; 
% params.sup_out = 0; %1 -> do not display iteration progress when running the algorithm. 
% times = zeros(2, length(8:20));
% count = 1; 
% for i = 8:20
%     disp('Iteration: ')
%     disp(i)
%     params.ns = i^2; %sources
%     
%     pre_processed_data = varproj_dataread(params); %extracts the data into a usable form
%     sub_data = varproj_sub(params, pre_processed_data); 
%     params = varproj_setparms(params, sub_data, pre_processed_data); 
%     
%     interpolated_cb = varproj_alg(params, sub_data, vec); 
%     %run FISTA (note options are for FISTA implementation)
% %     options.iterations = params.stop_crit*params.iter_crit*20; 
%     options.iterations = 2000; 
%     options.verbosity = 2; %all
%     interpolated_fista = fista_nucnorm_lapreg(sub_data,params,options);
%     
%     %times
%     times(1,count) = interpolated_cb.time;
%     times(2,count) = interpolated_fista.time; 
%     count=count+1; 
%     
% end

%% 2dinterp testing
% params.interp = 'res'; %'tt' or 'res'
% params.sub_used = 'ken'; %'ken' or 'rajiv' or 'real'
% %subsamp matrix you want
% params.norm_type = 'cb';%l2, l1, cb, lap note: lap still under construction
% pre_processed_data = varproj_dataread(params);
% params.switchdata = 'y';  %yes to switch data with model projected data, no to keep regular data and just get out stats
% [pre_processed_data, opstats] = varproj_optest(params, pre_processed_data); 
% %now do the subsampling step
% sub_data = varproj_sub(params, pre_processed_data); 
% params.sigma = 0; 
% params = varproj_setparms(params, sub_data, pre_processed_data); 
% params.sup_out = 0; %1 -> do not display iteration progress when running the algorithm. 
% interpolated = varproj_alg(params, sub_data, vec); 
% statsop = varproj_datavis(params, sub_data, interpolated, pre_processed_data); 

