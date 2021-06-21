function out = minFnc_run(params, data)

tic; 
%initialize algorithm parameters
k = params.k;  
ns = params.ns; 
Rx = params.Rx; 
Ry = params.Ry; 
interp = params.interp; 
lambda = params.lambda; 
g = params.g; 
 
%initialize data parameters
obs = data.obs; 
sub_tt = data.sub.tt; 
sub_res = data.sub.res; 
sub_u = data.sub.u;

%initialize
x = randn(Rx*sqrt(ns),k);
x = [x(:);x(:)]; 

Ao = opRestriction(Rx*Ry*ns,obs);%*sparse(diag(1./sub_u(:))); %option: don't include uncertainties
 
U = sparse(diag(1./sub_u(:))); 
%note that here we want to keep everything else in (RxRy)x(SxSy)
%space until the end since that's what the laplacian is created is

switch(interp)
   %pick what you are matching; also scale by reported uncertainties
   case{'tt'}
        b  = Ao*sub_tt(:);
   case{'res'}
        b  = Ao*sub_res(:); 
end


lbfgs_parms.Rx = Rx; 
lbfgs_parms.Ry = Ry; 
lbfgs_parms.ns = ns;
lbfgs_parms.b = b; 
lbfgs_parms.k = k; 
lbfgs_parms.gamma = g;
lbfgs_parms.lambda = lambda*100; 
lbfgs_parms.Lbig = params.Lbig;
lbfgs_parms.LtL = params.LaptLap; 
lbfgs_parms.AtA = Ao'*Ao;
lbfgs_parms.true = data.model.res; 

% mexAll; 
options = [];
options.display = 'iter';
options.maxFunEvals = 1000;
options.MaxIter = params.stop_crit*params.iter_crit; 
options.Method = 'lbfgs';
options.optTol = 1e-6; 
funobj = @(x)NLmatcomp(x, Ao, lbfgs_parms); 
[x, ~, ~, outputstruct] = minFunc(funobj,x,options);
e = Rx*sqrt(ns)*k;
L = x(1:e);
R = x(e+1:end);
L = reshape(L,Rx*sqrt(ns),k);
R = reshape(R,k, Ry*sqrt(ns));  
out.XM =  L*R;
out.time = toc; 
out.obj_val = outputstruct.trace.fval; 


end
