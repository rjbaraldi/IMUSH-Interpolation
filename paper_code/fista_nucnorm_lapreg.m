function info = fista_nucnorm_lapreg(data,params,options)
% FISTA_NUCNORM   Fast Iterative Shrinkage-Thresholding Algorithm
%
%    [X,INFO] = FISTA_NUCNORM(M,IDX,LAMBA,OPTIONS) solves
%
%       minimize 1/2 * ||X(idx)-M(idx)||_2^2 + lambda * ||X|| +
%       ||Lap(X)||_2^2
%
%    using the fast iterative shrinkage-thresholding algorithm
%    proposed by Beck and Teboulle. The following options are
%    provided:
%
%    OPTIONS
%       .iterations     maximum number of iterations
%       .verbosity      verbosity level: 0 = quiet, 1 = summary, 2 = all
%       .maxMatvec      maximum matrix-vector multiplies allowed
%       .maxRuntime     maximum runtime allowed (in seconds)
%       .prefix         prefix string
%

% Written by: Ewout van den Berg;
% September 28, 2009
%output param vector
LtL = params.LaptLap; 
k = params.k;  
ns = params.ns; 
Rx = params.Rx; 
Ry = params.Ry; 
interp = params.interp; 
lambda = params.lambda*50; 
g = params.g; 
I = params.I;
printevery = params.printevery;
suppress_out = params.sup_out; 
%initialize data parameters
obs = data.obs; 
sub_tt = data.sub.tt; 
sub_res = data.sub.res; 
sub_u = data.sub.u;

Ao = opRestriction(Rx*Ry*ns,obs)*I;%*sparse(diag(1./sub_u(:))); %option: don't include uncertainties
%preform matrix:
D = lambda*(Ao'*Ao)+g*LtL; 
U = spdiags(1./sub_u(:), 0, Rx*Ry*ns, Rx*Ry*ns); 
%note that here we want to keep everything else in (RxRy)x(SxSy)
%space until the end since that's what the laplacian is created is

switch(interp)
   %pick what you are matching; also scale by reported uncertainties
   case{'tt'}
        b  = Ao*sub_tt(:);
   case{'res'}
        b  = Ao*sub_res(:); 
end


% Start timer
t0 = tic;

% Set default parameters
iterations = 10000;
verbosity  = 1;
maxMatvec  = Inf;
maxRuntime = Inf;
prefix     = '';

if isfield(options,'iterations'), iterations = options.iterations; end
if isfield(options,'verbosity' ), verbosity  = options.verbosity;  end
if isfield(options,'maxMatvec' ), maxMatvec  = options.maxMatvec;  end
if isfield(options,'maxRuntime'), maxRuntime = options.maxRuntime; end
if isfield(options,'prefix'    ), prefix     = options.prefix;     end

% Estimate Lipschitz constant
L = 1/normest(params.Lbig)^2;
% Exit conditions (constants)
EXIT_ITERATIONS    = 1;
EXIT_OPTIMAL       = 2;
EXIT_MATVEC_LIMIT  = 3;
EXIT_RUNTIME_LIMIT = 4;


% ---------------------------------------------------------------------
% Log header
% ---------------------------------------------------------------------
if verbosity > 0
   fprintf('%s%s\n', prefix, repmat('-',40,1));
   fprintf('%s%5s  %9s  %9s  %9s\n',prefix,'iter','objective');
   fprintf('%s%s\n', prefix, repmat('-',40,1));
end

%Initialize
obj_value = []; 
X = randn(Rx*sqrt(ns),k);
X = X*X';  
d = lambda*Ao'*b; 
Xm = X;  t = 1; tm = t; 
iter = 1; status = 0;
tau = L; 
% =====================================================================
% MAIN LOOP
% =====================================================================
while 1
   Y = X + (tm - 1)/t*(X-Xm); %Y is a matrix here
   Xm = X; %update Xm
   tm = t; 
   [~,S,~] = svd(Y); %take svd
   s       = diag(S);   %get diagonals
   YNorm   = norm(s,1); %take norm of Y (nuke norm)
   G = Y(:) - tau^(-1)*(D*Y(:) - d); 
   f       = lambda*norm(Ao*X(:) - b,2)^2 + g*norm(params.Lbig*X(:),2)^2 + YNorm;
    feas = norm(Ao*X(:) - b,2); 
   % ------------------------------------------------------------------
   % Test exit conditions
   % ------------------------------------------------------------------
   if 2*iter >= maxMatvec
       status = EXIT_MATVEC_LIMIT;
   end
   if toc(t0) >= maxRuntime
       status = EXIT_RUNTIME_LIMIT;
   end
   if iter >= iterations
       status = EXIT_ITERATIONS;
   end
   

   % ------------------------------------------------------------------
   % Print log and act on exit conditions
   % ------------------------------------------------------------------
   if (mod(iter, printevery)==1) && ~suppress_out
            fprintf('iter: %d, ObjValue: %7.3f, feas: %7.3f\n', iter, f, feas);
   end

   if status ~= 0
      break;
   end
   
   % ------------------------------------------------------------------
   % Iterations begin here
   % ------------------------------------------------------------------
   iter   = iter + 1;
   [U,S,V]= svd(reshape(G, Rx*sqrt(ns), Ry*sqrt(ns)));
   sigma_g = max(diag(S) - 1/tau, 0);
   X      = U * diag(sigma_g) * V';
   t      = (1 + sqrt(1+4*t^2)) / 2;
   obj_value=[obj_value, norm(Y - data.model.res)]; 
end

% Determine status message
switch status
   case EXIT_ITERATIONS
      statusMsg = 'ERROR EXIT -- Too many iterations';
   case EXIT_OPTIMAL
      statusMsg = 'EXIT -- Optimal solution found';
   case EXIT_MATVEC_LIMIT
      statusMsg = 'EXIT -- Maximum matrix-vector operations reached';
   case EXIT_RUNTIME_LIMIT
      statusMsg = 'EXIT -- Maximum runtime reached';
   otherwise
      statusMsg = 'ERROR -- Unknown terminal condition';
end

if verbosity > 0
   fprintf('%s\n', prefix);
   fprintf('%s%s\n', prefix, statusMsg);
end

% Prepare output variables
info = struct();
info.XNorm     = YNorm; % [sic]
info.iter      = iter;
info.time = toc(t0);
info.status    = status;
info.statusMsg = statusMsg;
info.L         = L;
info.obj_val = obj_value; 
info.XM = Y;
