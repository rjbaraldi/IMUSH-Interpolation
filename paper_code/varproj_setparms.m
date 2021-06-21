function params = varproj_setparms(params,sub_data, pre_processed_data)
%% function used to set the parameters
%Input: takes in 1) preliminary parameter structure (used to make data)
%                2) sub_data -> subsampled data
%                3) pre_processed_data loaded data, used mostly for source
%                index
%Output: updated parameter structure. 
    

%%
%universal parameters
params.k = 40; 
params.converged = 1e-10; 
params.I = speye(params.Rx*params.Ry*params.ns); 
%case specific parameters
switch(params.norm_type)
    case{'l2', 'l1'}
        params.eta = .5;
        params.g = 0; %gamma (not used) in l1/l2
        params.printevery = 100; %iteration criteria
        params.stop_crit = 100;
        params.iter_crit = 5; 
        params.lambda = 1; 
%         params.etafact = 3; %how much you multiply by eta by after each iteration
        params.etafact = sum(svd(sub_data.real.res))/params.k; %roughly 4
    case{'cb','lap'}
        %note this syntax is different than in projection.pdf -> here,
        %gamma and eta are coefficients (not inverse coefficients, ie
        %1/(2*eta)
        params.eta = 1; 
        params.g = params.eta*(880)^2; %this goes in front of the laplacian
        params.nu = 1; 
        params.lambda = 45; 
        if params.sigma==0
            params.printevery= 10;
% %             stop_crit = 100;
%             params.stop_crit = 30; 
%             params.iter_crit = 14;
%             stop_crit = 100;
            params.stop_crit = 50; 
            params.iter_crit = 3;

        else
            params.printevery= 10;
            params.stop_crit = 30; 
            params.iter_crit = 3;
        end
        params.etafact = sum(svd(sub_data.real.res))/params.k;
        % Carl laplacian operator (remember this work on (rx*ry,sxsy) matrix which
        % is permuted compared to what we need for nuclear-norm A).
        [params.Lbig,~]  = makeL(params.Rx,params.Ry,pre_processed_data.indsort(1:params.ns)); %get laplacian (indsort refers to sources, so we can re-use preprocessed data
        temp = convertLap(1:size(params.Lbig,2), params.Rx, params.Ry, params.ns, 'index'); %get how index of row entries change
        P = opPermutation(temp(:)'); %make permutation operator based on how row entries change
        params.Lbig = sparse(params.nu*params.Lbig*P'); %apply permutation operator to laplacian to rearrange each row of the laplacian
        params.LaptLap = params.Lbig'*params.Lbig; %pre-form L'L
       
end



end