function out = varproj_sub(params, data)
%This function takes in some of the meta parameters that are being used to
%describe the grid, as well as the data in (RxRy)x(SxSy) form. Note that
%this function needs this form to do the subsampling correctly, since we
%are eliminating receiver positions across each source
sub_used = params.sub_used; 
Rx = params.Rx;
Ry = params.Ry;
perc = params.perc; 
ns = params.ns; 

real_tt = data.real.tt; %observed data
real_res= data.real.res; 
real_u = data.real.u; 
model_tt = data.model.tt; %model data
model_res = data.model.res; 
model_u = data.model.u; 

switch(sub_used)
    case{'ken', 'real'} 
        sub = find(real_tt(:)~=0); %find everything not zero
        [perc_samp, ~] = datasample(sub,floor(perc*numel(sub)), 'Replace', false); %sample without replacement
        sub_tt = real_tt(:);
        sub_tt(perc_samp) = 0; %replace indices of sampled data with zero
        sub_res = real_res(:); 
        sub_res(perc_samp)= 0; %repeat for all
        sub_u = real_u(:); 
        sub_u(perc_samp) = 0;
        sub_tt = reshape(sub_tt, Rx*Ry, ns); %reshape to get in RxRy x SxSy form
        sub_res = reshape(sub_res, Rx*Ry, ns); 
        sub_u = reshape(sub_u, Rx*Ry, ns);  
    case{'rajiv'}
        %in this case, we just sample across the RxRy x SxSy domain to get
        %rid of a station for all sources
        sub                  = randperm(Rx*Ry);
        sub                 = sub(1:floor(length(sub)*(1-perc)));
        sub_tt            = zeros(size(real_tt));
        sub_res            = zeros(size(real_tt));
        sub_u            = zeros(size(real_tt));
        sub_tt(sub,:)   = real_tt(sub,:);
        sub_res(sub,:)   = real_res(sub,:);
        sub_u(sub,:)   = real_u(sub,:);
        
end

%Sort columns via highest energy (can be tt or res, but if res use abs())
[idx,~] = max(abs(sub_res), [],1); %find the max entry in each column
[~,idx] = sort(idx, 2, 'descend'); %get indices for descending sort
%restructure the matrices via those column
sub_tt = sub_tt(:,idx); 
sub_res =sub_res(:,idx); 
sub_u = sub_u(:,idx); 
real_tt = real_tt(:,idx); 
real_res = real_res(:,idx); 
real_u = real_u(:,idx); 
model_tt = model_tt(:,idx); 
model_res = model_res(:,idx); 
model_u = model_u(:,idx); 


%reshape everything into RxSx x RySy plane
out.real.tt = reshape(permute(reshape(real_tt, Rx,Ry, sqrt(ns),sqrt(ns)),[1 3 2 4]),Rx*sqrt(ns), Ry*sqrt(ns)); 
out.real.res = reshape(permute(reshape(real_res, Rx,Ry, sqrt(ns),sqrt(ns)),[1 3 2 4]),Rx*sqrt(ns), Ry*sqrt(ns)); 
out.real.u = reshape(permute(reshape(real_u, Rx,Ry, sqrt(ns),sqrt(ns)),[1 3 2 4]),Rx*sqrt(ns), Ry*sqrt(ns)); 

out.sub.tt = reshape(permute(reshape(sub_tt, Rx,Ry, sqrt(ns),sqrt(ns)),[1 3 2 4]),Rx*sqrt(ns), Ry*sqrt(ns)); 
out.sub.res = reshape(permute(reshape(sub_res, Rx,Ry, sqrt(ns),sqrt(ns)),[1 3 2 4]),Rx*sqrt(ns), Ry*sqrt(ns)); 
out.sub.u = reshape(permute(reshape(sub_u, Rx,Ry, sqrt(ns),sqrt(ns)),[1 3 2 4]),Rx*sqrt(ns), Ry*sqrt(ns)); 

out.model.tt = reshape(permute(reshape(model_tt, Rx,Ry, sqrt(ns),sqrt(ns)),[1 3 2 4]),Rx*sqrt(ns), Ry*sqrt(ns)); 
out.model.res = reshape(permute(reshape(model_res, Rx,Ry, sqrt(ns),sqrt(ns)),[1 3 2 4]),Rx*sqrt(ns), Ry*sqrt(ns)); 
out.model.u = reshape(permute(reshape(model_u, Rx,Ry, sqrt(ns),sqrt(ns)),[1 3 2 4]),Rx*sqrt(ns), Ry*sqrt(ns)); 

out.no_obs = find(out.sub.tt(:)==0);%note these are indices, not binary representations
out.obs = find(out.sub.tt(:)~=0); 

%might also have to permute the variances from the dataread file

%plot tensors if need be
if params.figs_on
    varproj_figgen_field(real_res, 'Subsampled Residuals against 1D Model', 4, params, 'model_sub_rxrysxsy'); 
    varproj_figgen_field(out.model.res, 'Full Residuals against 1D Model', 2, params,'res_model_all'); 
    varproj_figgen_field(out.sub.res, 'Subsampled Residuals against 1D Model',4, params,'model_sub_rxsxrysy');
    varproj_figgen_field(out.model.res,'Single grid for 3D - 1D Model Residuals' , 3, params, 'res_model_single1') 

%     varproj_figgen_svd(real_res, 'Form 1', out.sub.res, 'Form 2', out.model.res, 'Model Truth')
%     figure;hold on
%     plot(svd(real_res)/max(svd(real_res)),'b.-');
%     plot(svd(out.sub.res)/max(svd(out.sub.res)),'r.-');
%     plot(svd(out.model.res)/max(svd(out.model.res)),'k.-');
%     ylabel('Relative Singular Value \sigma_i/\sigma_{max}', 'Fontsize', 22, 'Fontweight', 'Bold')
%     xlabel('ith index', 'Fontsize', 22, 'Fontweight', 'Bold')
%     legend('Form 1', 'Form 2', 'Model Truth'); 
%     title('Singular Value Decay', 'Fontsize', 22, 'Fontweight', 'Bold')
%     set(gca, 'Fontsize', 22, 'Fontweight', 'Bold')
%     hold off
end
end
