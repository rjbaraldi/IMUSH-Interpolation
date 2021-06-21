function stats = varproj_datavis(params,subdata, interpolated, real_data, normtitle, savekey)
%Input: completed data from algorithm

%Output: Figures for residual and travel times, as well as a latex table
%row of parameters used and performance metrics. 


%initialize alg parameters
sub_used = params.sub_used; 
norm_type = params.norm_type; 
interp = params.interp; 
sigma = params.sigma; 
perc = params.perc; 
g = params.g; 
eta = params.eta;
ns = params.ns; 
Rx = params.Rx; 
Ry = params.Ry; 
model_cutoff = params.model_cutoff; 

%initialize data parameters
%subsampled indices
obs = subdata.obs; 
no_obs = subdata.no_obs;
%objective value
obj_value = interpolated.obj_val; 
%subsampled data
sub_tt = subdata.sub.tt; 
sub_res = subdata.sub.res;
%model data
model_tt = subdata.model.tt; 
model_res = subdata.model.res;
%real data
%note that for model/ken - real data is noisy model data and if 
%perc=0, then subdata.real.res = subdata.sub.res
real_res = subdata.real.res; 
real_tt = subdata.real.tt; 
%uncertainties
sub_u = subdata.sub.u;
XM = interpolated.XM; 
timetook = interpolated.time;  

if params.figs_on
switch(interp)
    case{'tt'} %get travel time results
        rec_tt = XM; 
        rec_resid = sub_tt-XM;
        rec_resid(no_obs)=0;
        %now consider travel times
        varproj_figgen_field(sub_tt,strcat('True/observed TTs: ',normtitle),2, params, strcat('obs_',savekey,'_all'));
        varproj_figgen_field(rec_tt,strcat('Interpolated TTs: ',normtitle),2, params,strcat('interpolated_',savekey,'_all'));
        varproj_figgen_field(rec_resid,strcat('Difference: ',normtitle),2, params,strcat('difference_',savekey,'_all')); 
%         varproj_figgen_field(model_tt,strcat('Model TTs: ',normtitle),2, params,strcat('obs_',savetitle,'_all')); 
        varproj_figgen_field(sub_res, strcat('True/observed Res: ',normtitle),2, params,strcat('obs_',savekey,'_all')); 
        varproj_figgen_field(rec_resid,strcat('Interpolated Res: ',normtitle),2, params,strcat('interpolated_',savekey,'_all')); 
        varproj_figgen_field(sub_tt,strcat('True/observed TTs: ',normtitle),3, params,strcat('obs_',savekey,'_single1'));
        varproj_figgen_field(rec_tt,strcat('Interpolated TTs: ',normtitle),3, params,strcat('interpolated_',savekey,'_single1'));
        varproj_figgen_field(rec_resid,strcat('Difference: ',normtitle),3, params,strcat('difference_',savekey,'_single1')); 
%         varproj_figgen_field(model_tt,strcat('Model TTs: ',normtitle),3, params,strcat('obs_',savekey,'_all')); 
        varproj_figgen_field(sub_res, strcat('True/observed Res: ',normtitle),3, params,strcat('obs_',savekey,'_single1')); 
        varproj_figgen_field(rec_resid,strcat('Interpolated Res: ',normtitle),3, params,strcat('obs_',savekey,'_single1'));
        figure; %plot difference
        if perc==0
            varproj_figgen_field(sub_res-rec_resid,strcat('Difference in Observed: ',normtitle),2, params,strcat('difference_',savekey,'_all')); 
            varproj_figgen_field(sub_res-rec_resid,strcat('Difference in Observed: ',normtitle),3, params,strcat('difference_',savekey,'_single1')); 
        else
            tr = zeros(size(real_tt));
            tr(sub_res==0 & real_res~=0) = real_res(sub_res==0 & real_res~=0); 
            trec = zeros(size(rec_resid)); 
            trec(sub_res==0 & real_res~=0) = rec_resid(sub_res==0 & real_res~=0); 
            varproj_figgen_field(tr - trec,strcat('Difference in Recovered: ', normtitle),2, params,strcat('difference_',savekey,'_all')); 
            varproj_figgen_field(tr - trec,strcat('Difference in Recovered: ', normtitle),3, params,strcat('difference_',savekey,'_single1')); 
        end
        
       
%         varproj_figgen_field(model_res, 'Full Model Res',2, params,strcat('model_',savekey,'_all'));
%         varproj_figgen_field(model_res, 'Full Model Res',3, params,strcat('model_',savekey,'_single1'));
        %plot SVD info
%         varproj_figgen_svd(XM, 'Interpolated TT', sub_tt, 'Observed TT', model_tt, 'Model TT')

    case{'res'}
        rec_resid = XM; 
        rec_resid(no_obs)=0;
        temp = sub_res;
        rec_tt = XM - (sub_res - sub_tt); 
        rec_tt(no_obs)=0;
        
        varproj_figgen_field(sub_res,strcat('Observed Res: ', normtitle),4, params,strcat('obs_',savekey,'_all')); 
        varproj_figgen_field(XM,strcat('Interpolated Res: ', normtitle),2, params,strcat('interpolated_',savekey,'_all')); 
        varproj_figgen_field(sub_res,strcat('Observed Res in Single Source: ', normtitle),3, params,strcat('obs_',savekey,'_single1')); 
        varproj_figgen_field(XM,strcat('Interpolated: ', normtitle),3, params,strcat('interpolated_',savekey,'_single1')); 
        figure; %plot difference
%         if perc==0
        if strcmp(sub_used, 'ken')
            varproj_figgen_field(abs(model_res-XM),strcat('True Difference: ', normtitle),2, params, strcat('difference_',savekey,'_all')); 
            varproj_figgen_field(abs(model_res-XM),strcat('True Difference in Single Source: ',normtitle),3, params, strcat('difference_',savekey,'_single1')); 
        
        else
            tr = zeros(size(real_tt));
            tr(sub_res==0 & real_res~=0) = real_res(sub_res==0 & real_res~=0); 
            trec = zeros(size(rec_resid)); 
            trec(sub_res==0 & real_res~=0) = rec_resid(sub_res==0 & real_res~=0); 
            varproj_figgen_field(tr - trec,strcat('Recovered Difference: ', normtitle), 2, params, strcat('difference_',savekey,'_all')); 
            varproj_figgen_field(tr - trec,strcat('Recovered Difference (single): ',normtitle), 3, params, strcat('difference_',savekey,'_single1')); 
        end
%         varproj_figgen_field(model_res,'Model Res',2, params); 
%         varproj_figgen_field(model_res,'Model Res in Single Source',3, params); 
        
        %now consider travel times
%         varproj_figgen_field(sub_tt, 'True/observed TTs',2, params); 
%         varproj_figgen_field(rec_tt, 'Interpolated TTs',2, params);
%         
%         varproj_figgen_field(temp-rec_resid,'Difference', 2, params); 
%         varproj_figgen_field(model_tt,'Model TTs',2, params); 
%         
        
%         varproj_figgen_svd(XM, 'Interpolated Res', sub_res, 'Observed Res', model_res, 'Model Res')
        
        
        
end
%plot objective function decay
% figure;
% plot(1:numel(obj_value), obj_value, 'k', 'LineWidth', 2)
% xlabel('Iteration', 'Fontsize', 22, 'Fontweight', 'Bold')
% ylabel('Objective Function Value', 'Fontsize', 22, 'Fontweight', 'Bold')
% title('Objective Function Decay', 'Fontsize', 22, 'Fontweight', 'Bold')
% set(gca, 'Fontsize', 22, 'Fontweight', 'Bold')
end
switch(interp)
    case{'tt'}       
        rec_tt = XM; 
        rec_resid = sub_tt-XM;
    case{'res'}
        rec_resid = XM; 
        rec_tt = zeros(size(XM)); 
        rec_tt(real_tt~=0) = XM(real_tt~=0) - (-real_tt(real_tt~=0)+ real_res(real_res~=0));      
end
%other metrics to use
switch(sub_used)
    case{'rajiv', 'ken'}
        stats.true_recov = norm(model_res(obs) - rec_resid(obs))/sqrt(numel(obs));%same as rms_seen
        stats.unseen_recov = norm(model_res(no_obs) - rec_resid(no_obs))/sqrt(numel(no_obs));%same as rms_unseen for
%         the whole grid
        stats.SNRtt = -20*log10(norm(model_tt(obs)-rec_tt(obs),'fro')/norm(real_tt(obs),'fro'));
        stats.SNRre = -20*log10(norm(model_res(obs)-rec_resid(obs),'fro')/norm(real_res(obs),'fro'));

    case{'real'}
        op = real_data.op;
        load('EV_135_syn6.mat');
        EV_used = EV(real_data.indsort(1:ns));
        stats.meanrec = zeros(1,ns); 
        stats.medianrec = zeros(1, ns); 
        stats.true_recov_data = zeros(1,ns);
        stats.unseen_recov = zeros(1, ns); 
        stats.num_obs = zeros(1, ns);
        stats.num_noobs = zeros(1, ns);
        stats.true_recov_data = sqrt(sum(stats.true_recov_data))/sum(stats.num_obs); 
        stats.true_recov = norm(real_res(obs) - rec_resid(obs))/sqrt(numel(obs));%same as rms_seen
        stats.unseen_recov = norm(real_res(no_obs) - rec_resid(no_obs))/sqrt(numel(no_obs));%same as rms_unseen
        stats.SNRtt = -20*log10(norm(real_tt(real_tt~=0)-rec_tt(real_tt~=0),'fro')/norm(rec_tt(real_tt~=0),'fro'));
        stats.SNRre = -20*log10(norm(real_res(real_res~=0)-rec_resid(real_res~=0),'fro')/norm(rec_resid(real_res~=0),'fro'));
        rec_resid = reshape(permute(reshape(rec_resid, Rx,Ry, sqrt(ns),sqrt(ns)),[1 3 2 4]),Rx*Ry,ns);
        stats.disttocent=[]; %center is [95, 95]
        stats.cumulabsdiff=[]; 
        for i = 1:params.ns 
            resnew = rec_resid(:,i)'*op{i};
            resold = EV_used(i).resid(model_cutoff:end)'; 
            for j = 0:numel(EV_used(i).resid)-model_cutoff
                stats.disttocent=[stats.disttocent ...
                    sqrt((EV_used(i).x(model_cutoff+j)-95)^2+((EV_used(i).y(model_cutoff+j)-95)^2))]; 
                stats.cumulabsdiff=[stats.cumulabsdiff abs(resnew(j+1)-resold(j+1))];
            end
            stats.num_obs(i) = size(model_cutoff:numel(EV_used(i).resid),2);
            stats.rmse(i) = sum((resnew - resold).^2)/stats.num_obs(i); 
            stats.meanrec(i) = mean(abs(resnew - resold)); 
            stats.medianrec(i) = median(abs(resnew -resold)); 
        end
end
%print all this stuff out for 
fprintf('Experiment: %s/%s/%s & %.2f & %.0f & %.2f & %.2f & %.1e & %.1e & %.2f & %.2f & %.3f & %.3f \n', norm_type,interp,sub_used,sigma,numel(find(sub_res~=0))/numel(real_res)*100,...
    timetook, g, eta, stats.SNRtt, stats.SNRre, stats.true_recov, stats.unseen_recov); 
% size(find(sub_tt(:)~=0),1)/numel(real_tt)*100,
end