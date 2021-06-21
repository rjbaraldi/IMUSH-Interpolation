function [results, interpgrids, bootstats, finalstats] = varproj_bootstrapping(params, data, vec)
%This function does the bootstrapping for random sampling percentage
%specified by params.perc. For now (3/4/18), bootstrapping data visualization is
%done here instead of varproj_datavis due to the potential size of the
%data. Also, it is recommended that this is not run on a personal computer
%since the size of the data stored can blow up quickly. This may be changed
%in the future; 

%Input: params-> algorithm parameter structure
%       data -> pre_processed data structure; note this has to be the raw
%       data read since the routine varproj_sub() does the subsampling
%       stuff. 
%       vec -> vector function
%Output: bootstruct-> data structure that contains the data interpolated by
%the bootstrapping runs, as well as the mean interpolated tensor and the
%variance of the interpolated tensor. 


%Initialize algorithm parameters
Rx = params.Rx; 
Ry = params.Ry; 
ns = params.ns; 
interp = params.interp;
sub_used = params.sub_used; 
bootruns = params.bootruns; 
figs_on = params.figs_on;
save_on = params.save_on; 
printevery = params.printevery; 
%data structure for storing results
results = cell(1,bootruns);
bootstats = zeros(4, bootruns); %1 is true recovered, 2 is unseen, 3 is snrt, 4 is snrr
%tensor for doing statistics
interpgrids = zeros(Rx*sqrt(ns), Ry*sqrt(ns), bootruns); 


parfor i = 1:bootruns
    %Repeat the subsampling step with the same percentage removed
    temp_data = varproj_sub(params, data); 
    results{i} = varproj_alg(params, temp_data, vec);
    interpgrids(:,:,i) = results{i}.XM; 
    if(mod(i, printevery)==1)
            fprintf('Bootstrap Iter: %d/%d\n', i, bootruns);
    end
    
    %also compute performance metrics on mean/var
    switch(interp)
        case{'tt'}       
            rec_tt = results{i}.XM; 
            rec_resid = temp_data.sub.tt-results{i}.XM;
        case{'res'}
            rec_resid = results{i}.XM; 
            rec_tt = results{i}.XM - (temp_data.real.res - temp_data.real.tt);      
    end
    %define parameters for statistics
    obs = temp_data.obs; 
    no_obs = temp_data.no_obs; 
    real_res = temp_data.real.res;
    real_tt = temp_data.real.tt; 
    model_res = temp_data.model.res;
    sub_res = temp_data.sub.res; 
    sub_tt = temp_data.sub.tt; 
    %other metrics to use
    switch(sub_used)
        case{'rajiv', 'ken'}
        bootstats(:,i) =[norm(real_res(temp_data.obs) - rec_resid(obs))/sqrt(numel(obs));...
            norm(model_res(no_obs) - rec_resid(no_obs))/sqrt(numel(no_obs));...
            -20*log10(norm(real_tt(obs)-rec_tt(obs),'fro')/norm(real_tt(obs),'fro'));...
            -20*log10(norm(real_res(obs)-rec_resid(obs),'fro')/norm(real_res(obs),'fro'))];

        case{'real'}
           bootstats(:,i) =[norm(real_res(temp_data.obs) - rec_resid(obs))/sqrt(numel(obs));...
               norm(real_res(sub_res==0 & real_res~=0) - rec_resid(sub_res==0 & real_res~=0))/...
                sqrt(numel(results{i}.XM(sub_res==0 & real_res~=0)));...
                -20*log10(norm(real_tt(sub_tt==0 & real_tt~=0)-rec_tt(sub_tt==0 & real_tt~=0),'fro')/norm(rec_tt(sub_tt==0 & real_tt~=0),'fro'));...
                -20*log10(norm(real_res(sub_res==0 & real_res~=0)-rec_resid(sub_res==0 & real_res~=0),'fro')/norm(rec_resid(sub_res==0 & real_res~=0),'fro'))];

    end

    
    
    
    
end

finalstats.mean = mean(interpgrids, 3); 
finalstats.var = var(interpgrids,[], 3); 
if figs_on
    %plot mean
    figure; 
    imagesc(finalstats.mean); 
    xlabel('Mean - Interpolated Resids', 'Fontsize', 22, 'Fontweight', 'Bold'); colorbar;
    set(gca, 'Fontsize', 22, 'Fontweight', 'Bold');
    
    %plot variance
    figure; 
    imagesc(finalstats.var); 
    xlabel('Variance - Interpolated Resids', 'Fontsize', 22, 'Fontweight', 'Bold'); colorbar;
    set(gca, 'Fontsize', 22, 'Fontweight', 'Bold'); 
end


%save if you want
if save_on
   savestr = ['bootdata','_',date]; 
   save(savestr, 'results', 'bootstats', 'finalstats'); 
end

end

