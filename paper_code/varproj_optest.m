function [out, stats] = varproj_optest(params,real_data)
%initialize everything
ns = params.ns;  
model_cutoff = params.model_cutoff; 
Rx = params.Rx;
Ry = params.Ry;
Gx = params.Gx; 
Gy = params.Gy; 
% perc = params.perc; 
resolution = params.resolution; 
load('EV_135_syn6.mat');
EV_used = EV(real_data.indsort(1:ns));
stats.meanrec = zeros(1,ns); 
stats.medianrec = zeros(1, ns); 
stats.true_recov_data = zeros(1,ns);
stats.unseen_recov = zeros(1, ns); 
stats.num_obs = zeros(1, ns);

stats.proj_rms_data = zeros(1, ns);
stats.proj_eff_data = cell(1, ns);
stats.proj_rms_model = zeros(1,ns); 
stats.proj_eff_model = cell(1,ns); 
if strcmp(params.switchdata, 'y')
    out.indsort = real_data.indsort; 
    out.model.tt = real_data.model.tt; 
    out.model.res = real_data.model.res; 
    out.model.u = real_data.model.u;
    real_tt = zeros(Rx*Ry,ns);
    real_res = zeros(Rx*Ry, ns);
    real_u = zeros(Rx*Ry,ns);
    out.op = cell(1,ns); 
end
for i = 1:ns
    %1) we want to evaluate the efficacy of the operator, so we can first
    %project from real data to a grid and back again, measuring the error
    coord = [EV_used(i).x(model_cutoff:end)', EV_used(i).y(model_cutoff:end)']; %coordinates used if you need them
    op = opLInterp2D(Gx(1):resolution:Gx(2), Gy(1):resolution:Gy(2), coord)';%*eye(size(coord,1));
    %project onto and off of grid again, then subtract
    %turn op into a matrix
    O = zeros(size(op)); 
    for kk = 1:size(O,2)
        t = zeros(size(op,2),1); 
        t(kk) = 1; 
        O = op*t; 
    end
%     stats.proj_eff_data{i} = op'*(real_data.real.res(:,i)) - EV_used(i).resid(model_cutoff:end);
    stats.proj_eff_data{i} = (O'*O)\(O'*real_data.real.res(:,i)) - EV_used(i).resid(model_cutoff:end); 
    stats.proj_rms_data(i) = norm(stats.proj_eff_data{i})/numel(stats.proj_eff_data{i}); 
    
    %do the same thing for the model
    x = op*(op'*(double(real_data.model.res(:,i)))); 
    stats.proj_eff_model{i} = x - real_data.model.res(:,i); 
    stats.proj_rms_data(i) = norm(stats.proj_eff_model{i})/numel(stats.proj_eff_model{i}); 
    if strcmp(params.switchdata, 'y')
       out.op{i} = op;  
       x = op*(op'*(double(real_data.model.tt(:,i))));
       var_temp = double(diag(1./EV_used(i).uncert(model_cutoff:end)));
       prod1 = zeros(model_cutoff-1, size(var_temp,1)); 
       for k = 1:size(var_temp,1)
           prod1(:,k) = op*var_temp(:,k); 
       end
       out.var{i} = prod1*op';         
       real_tt(:,i)= x;
       real_res(:,i)= op*(op'*(double(real_data.model.res(:,i))));
       real_u(:,i)= op*(op'*(double(real_data.model.u(:,i))));       
       SEM = std(x)/sqrt(length(x));               % Standard Error
       ts = tinv([0.01  0.99],length(x)-1);      % T-Score
       CI = mean(x) + ts*SEM;
       int_noise = 2*CI(2); 
%            int_noise = min(EV_used(i).tt(model_cutoff:end)); 
            %get rid of noise
       opidx = real_tt(:,i)<int_noise & real_tt(:,i)>-int_noise;
       real_tt(opidx, i)=0;
       real_u(opidx,i)=0;
       real_res(opidx,i)=0; 
    end
end
if strcmp(params.switchdata, 'y')
    out.real.tt = real_tt; 
    out.real.res = real_res; 
    out.real.u = real_u; 
else
    out =real_data; 
end



end
    
    
    
    
    
    



