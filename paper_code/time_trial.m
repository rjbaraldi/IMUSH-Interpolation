%% Preliminary data sorting - again, you should only have to run the last section
%beware the last section takes a long time to build the laplacian for each.
%
clear; clc; %close all
%define vec function
vec = @(x) x(:);
%add paths you need
addpath(genpath('../pSPOT'));
addpath(genpath('../spot'));
addpath(genpath('./minFunc_2012')); 

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
params.model_cutoff = 400+1; %build in buffer -> only used when manually binning data.
params.perc = .2; %note: perc = 0 => see everything. perc = 1 => see nothing
params.figs_on = 0; %1 -> figs on, 0 -> figs off
params.interp = 'res'; %'tt' or 'res'
params.sub_used = 'ken'; %'ken' or 'rajiv' or 'real'
%subsamp matrix you want
params.norm_type = 'cb';%l2, l1, cb, lap note: lap still under construction
pre_processed_data = varproj_dataread(params); %extracts the data into a usable form
sub_data = varproj_sub(params, pre_processed_data);  
params.sigma = 0; 
params = varproj_setparms(params, sub_data, pre_processed_data);
%% FISTA and variable relaxation comparison
%-doesn't work
params.sup_out = 0; %1 -> do not display iteration progress when running the algorithm. 
times = zeros(2, length(8:20));%ie mat comp sizes 64 - 400
count = 1; 
for i = 8:20
    disp('Iteration: ')
    disp(i)
    params.ns = i^2; %sources
    
    pre_processed_data = varproj_dataread(params); %extracts the data into a usable form
    sub_data = varproj_sub(params, pre_processed_data); 
    params = varproj_setparms(params, sub_data, pre_processed_data); 
    
    interpolated_cb = varproj_alg(params, sub_data, vec); 
    %run FISTA (note options are for FISTA implementation)
%     options.iterations = params.stop_crit*params.iter_crit*20; 
    options.iterations = 2000; 
    options.verbosity = 2; %all
    interpolated_fista = fista_nucnorm_lapreg(sub_data,params,options);
    
    %times
    times(1,count) = interpolated_cb.time;
    times(2,count) = interpolated_fista.time; 
    count=count+1; 
    
end


%% SKIP THIS SECTION FOR NOW
% nsrcx = 1; 
% nsrcy = nsrcx; 
% nrecx = 20;
% nrecy = 20;
% f = @(x,y) exp(x.^2+y.^2).*sin(pi.*x.*y);
% x = linspace(0,1,nrecx);
% y = linspace(0,1,nrecy); 
% [X, Y] = meshgrid(x,y);
% z = f(X,Y); 
% % clear out; 
% out = reshape(repmat(z,[1, nsrcx*nsrcy]), nrecx, nrecy, nsrcx*nsrcy);
% true = out; 
% mode.noise = 1;
% if mode.noise==1
%     tot_noise = 0; 
%     perc = .8; %used to be 0.01
%     idx = setdiff(1:nsrcx*nsrcy,indm);%set of guy's you've observed
%     indn = randperm(length(idx)); %random permutation of them
%     r = floor(length(indn)*perc); %define length for ease
%     %here, we corrupt with really strong points of noise  
%     for i = 1:r
%         switch(dataclass)
%             case{'l2'}
%                 alpha = 1e-2;
%                 normb = norm(out(:,:,idx(indn(i))),2);
%                 noisemat = alpha*normb*(randn(nrecx, nrecy));
%             case{'l1','l0'}
%                 %here we want to add very large noise values to a few
%                 %points
%                 alpha = 1e-1; 
%                 normb = norm(out(:,:,idx(indn(i)),1)); 
%                 noisemat = vec(alpha*normb*(randn(nrecx, nrecy))); %formerly 1e3
%                 [~, nind] = sort(abs(noisemat), 'descend');
%                 noisemat(nind(2+1:end))=0; %keep everything but the strongest
%                 noisemat=reshape(noisemat, nrecx, nrecy);
%             case{'linf'}
%                 alpha = 1e-3; 
%                 normb = norm(out(:,:,idx(indn(i))),inf); 
%                 noisemat = real(alpha*normb*exp(randn([nrecx, nrecy])*1i));
%         end
%         out(:,:,idx(indn(i))) = out(:,:,idx(indn(i))) + noisemat; %add back to values? 
%     end
% end
% out = reshape(permute(reshape(out,nrecx,nrecy,nsrcx,nsrcy),[1 3 2 4]),nrecx*nsrcx,nrecy*nsrcy);
% true = reshape(permute(reshape(true,nrecx,nrecy,nsrcx,nsrcy),[1 3 2 4]),nrecx*nsrcx,nrecy*nsrcy);

% params.ns = nsrcy*nsrcx; 
% params.Rx = nrecx; 
% params.Ry = nrecy; 
% params.interp = 'res'; %'tt' or 'res'
% params.sub_used = 'ken'; %'ken' or 'rajiv' or 'real'
% %subsamp matrix you want
% params.norm_type = 'cb';%l2, l1, cb, lap note: lap still under construction
% params.sigma = 0; 
% params = varproj_setparms(params, out); 
% params.sup_out = 0; %1 -> do not display iteration progress when running the algorithm. 
% b = out;
% ind = find(out==0);
% % k = rank(z);
% k = 3; 
% params.k = k; 
% params.eta = .1; 
% count = 0; 
% params.printevery = 50;
% stop_crit = 200; 
% params.stop_crit = stop_crit; 
% params.iter_crit = 3;
% params.method = 'nondist';
% params.obs = find(b(:)~=0); 
% params.no_obs = find(b(:)==0);
% params.converged = 1e-10;
% params.numr       = nsrcx*nrecx;
% params.numc       = nrecy*nsrcy;
% params.nr         = k;
% params.ind        = ind;
% params.mode       = 1;
% params.ls         = 1;
% % t=sum(svd(out))/params.k; 
% params.L = 1*randn(nsrcx*nrecx,k);
% params.R = 1*randn(nrecy*nsrcy,k);
% params.Ao = opRestriction(nsrcx*nrecx*nsrcy*nrecy, params.obs);
% params.An = opRestriction(nsrcx*nrecx*nsrcy*nrecy, params.no_obs);
% 
% svdtrial = 10; 
% times = zeros(2, svdtrial); 
% for i = 1:svdtrial
%     interpolated_cb = varproj_alg(params, sub_data, vec); 
%     %run FISTA (note options are for FISTA implementation)
%     options.iterations = params.stop_crit*params.iter_crit*10; 
%     options.verbosity = 2; %all
%     interpolated_fista = fista_nucnorm_lapreg(sub_data,params,options);
%     times(1, i) = interpolated_cb.time; 
%     times(2, i) = interpolated_fista.time; 
% end
% 
% 
