function out = varproj_dataread(params, dataset)
%%%%%
%code for reading/parsing the data
%Inputs: in -> struct that contains Rx, Ry, Gx, Gy, resolution, ns
%        dataset -> dataset you want to load
%Outputs: matrices you'll perform interpolation on for the model, data, and
%subsampled.
ns = params.ns; 
sub_used = params.sub_used; 
model_cutoff = params.model_cutoff; 
Rx = params.Rx;
Ry = params.Ry;
Gx = params.Gx; 
Gy = params.Gy; 
% perc = params.perc; 
resolution = params.resolution; 
if nargin<2
    load('EV_135_syn6.mat');
else
    %make sure dataset has same syntax
    out = load(dataset);
    EV = out.EV; 
    evind = out.evind; 
end

%% Initialize stuff
%get highest number of rec/source
num_rec_per_src = zeros(1,numel(EV));
evx = zeros(1,numel(EV));
evy = zeros(1,numel(EV)); 
for i = 1:numel(EV) %for # of events
    evx(i) = EV(i).evx; %get x for source
    evy(i) = EV(i).evy; %get y for source
    num_rec_per_src(i) = numel(EV(i).y); %get number of recievers per source
end
[~,out.indsort] = sort(num_rec_per_src,2, 'descend'); % use this sorting to load the data, this sorting make sure that
% you read the events in correct order according to the spatial grid (old)
%pick out top 1-ns number of sources
EV_used = EV(out.indsort(1:ns)); 
if params.figs_on==1
    temp=Gx(1):resolution:Gx(2); 
    x = EV_used(1).x(model_cutoff:end);
    y = EV_used(1).y(model_cutoff:end);
    x1 = x(Gx(1)<x & Gx(2)>x & y>Gx(1) & y<Gx(2)); 
    y1 = y(Gx(1)<x & Gx(2)>x & y>Gx(1) & y<Gx(2)); 
    figure;
    hold on;
    plot(evx(out.indsort(1:ns)), evy(out.indsort(1:ns)), '*b', 117.5, 117.5, '^r',evx(out.indsort(1)), evy(out.indsort(1)),'*r', 'LineWidth', 2);
    plot(temp, Gx(1)*ones(size(temp)), '--k')
    plot(Gx(1)*ones(size(temp)), temp, '--k')
    plot(Gx(2)*ones(size(temp)), temp, '--k')
    plot(temp, Gx(2)*ones(size(temp)), '--k')
    xlabel('Source - x (km)', 'Fontsize', 22, 'fontweight', 'bold'); 
    ylabel('Source - y (km)', 'fontsize', 22, 'fontweight', 'bold'); 
    title('Source Positions around Mt. St. Helens', 'Fontsize', 22, 'fontweight', 'bold');
    set(gca, 'fontsize', 22, 'fontweight', 'bold');
    axis([0 215 0 215])
    daspect([1 1 1])
    hold off
    [X, Y] = meshgrid(temp, temp);
    figure;
    hold on;
    plot(117.5, 117.5, '^r', 'LineWidth', 2);
    plot(X, Y, '*k')
    plot(x1, y1, '.b', 'MarkerSize', 10)
    xlabel('R_x (km)', 'Fontsize', 22, 'fontweight', 'bold'); 
    ylabel('R_y (km)', 'fontsize', 22, 'fontweight', 'bold'); 
    title('Station Positions around Mt. St. Helens', 'Fontsize', 22, 'fontweight', 'bold');
    set(gca, 'fontsize', 22, 'fontweight', 'bold');
    grid on
    axis([65 170 65 170])
    daspect([1 1 1])
    hold off
end




model_tt = zeros(Rx*Ry, ns); 
model_res = zeros(Rx*Ry,ns); 
model_u = zeros(Rx*Ry,ns);
real_tt = zeros(Rx*Ry,ns);
real_res = zeros(Rx*Ry, ns);
real_u = zeros(Rx*Ry,ns);

switch sub_used
    %in the rajiv case, we are sampling a percentage of the full grid
    %So here, we are using the model data as 'true' data
    case 'rajiv'
        for i = 1:ns
            model_tt(:,i)= EV_used(i).tt(1:model_cutoff-1);
            model_res(:,i)= EV_used(i).resid(1:model_cutoff-1);
            model_u(:,i)= EV_used(i).uncert(1:model_cutoff-1);
        end
        real_tt = model_tt;
        real_res = model_res;
        real_u = model_u; 
    case 'ken'
        %in the 'ken' case, we are using model data at the positions close
        %to the receivers. the subsampling is around 15% for each source
        %Carl has provided custom indices for this
        %we also take the model data for compariso 
        for i = 1:ns
            ins                = EV_used(i).indsub;
            real_tt(ins,i) = EV_used(i).ttsub;
            real_res(ins,i)  = EV_used(i).residsub;
            real_u(ins,i) = EV_used(i).uncertsub; 
            model_tt(:,i)  = EV_used(i).tt(1:model_cutoff-1);
            model_res(:,i)     = EV_used(i).resid(1:model_cutoff-1);
            model_u(:,i)  = EV_used(i).uncert(1:model_cutoff-1);
        end
        
    case 'real'
        %In this, we use the true data exactly and put them into gridpoints
        %assigned by receiver coordinates
        %the commented out code is used to the the same for the model, but
        %note that in the current grid and resolution, the 1:400 model selection is
        
        %Note here that we are no longer subsampling receiver positions
        %from off the grid (ie all receiver positions are taken into
        %account)
        out.op = cell(1, ns); 
        out.var = cell(1,ns); 
        for i = 1:ns
            coord = [EV_used(i).x(model_cutoff:end)', EV_used(i).y(model_cutoff:end)'];
            Op = opLInterp2D(Gx(1):resolution:Gx(2), Gy(1):resolution:Gy(2), coord)';
            %this creates a lot of extra points that we would then
            %subsample, so we have to trim the fat
            out.op{i} = Op; 
            x = Op*double(EV_used(i).tt(model_cutoff:end)');
            %random work-around to escape error in sparse applying data
            var_temp = double(diag(1./EV_used(i).uncert(model_cutoff:end)));
            prod1 = zeros(model_cutoff-1, size(var_temp,1)); 
            for k = 1:size(var_temp,1)
                prod1(:,k) = Op*var_temp(:,k); 
            end
            out.var{i} = prod1*Op'; 
           
            real_tt(:,i)= x;
            real_res(:,i)= Op*double(EV_used(i).resid(model_cutoff:end));
            real_u(:,i)= Op*double(EV_used(i).uncert(model_cutoff:end));
            
            SEM = std(x)/sqrt(length(x));               % Standard Error
            ts = tinv([0.01  0.99],length(x)-1);      % T-Score
            CI = mean(x) + ts*SEM;
            int_noise = 2*CI(2); 
%             int_noise = min(EV_used(i).tt(model_cutoff:end)); 
            %get rid of noise
            opidx = real_tt(:,i)<int_noise & real_tt(:,i)>-int_noise;
            real_tt(opidx, i)=0;
            real_u(opidx,i)=0;
            real_res(opidx,i)=0; 
            %if current grid is [80 180]/5, then can just use this for the
            %model
            model_tt(:,i)  = EV_used(i).tt(1:model_cutoff-1);
            model_res(:,i)     = EV_used(i).resid(1:model_cutoff-1);
            model_u(:,i)  = EV_used(i).uncert(1:model_cutoff-1);
         end
end
% 
% if params.figs_on
%     ax = Gx(1):resolution:Gx(2);
%     %if need be, plot subsampled data
%     figure; imagesc(ax, ax, flipud(imrotate(reshape(model_res(:,1),Rx, Ry), 90))); colorbar; xlabel('Receiver - x (km)', 'Fontsize', 22, 'Fontweight', 'Bold')
%     ylabel('Receiver - y (km)', 'Fontsize', 22, 'Fontweight', 'Bold')
%     title('Single Source-Receiver grid for 1d - 3d Model Residuals', 'Fontsize', 22, 'Fontweight', 'Bold')
%     set(gca, 'Fontsize', 22, 'Fontweight', 'Bold','YDir','normal');
% end

%store out structure
out.model.tt = model_tt; 
out.model.res= model_res; 
out.model.u = model_u; 
out.real.tt = real_tt; 
out.real.res = real_res; 
out.real.u = real_u; 



end