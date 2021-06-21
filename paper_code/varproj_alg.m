%% Initialize algorithm by cases
%Initialize w
function out = varproj_alg(params, data, vec)
%%%variable projection algorithm
%Input: 
%       params -> structure containing the parameters used in the low-rank min. Note
%           that some of the combined formulation also has pre-formed matrices.
%           Note that LaptLap is formed in the wrapper in passed in, but the restriction
%           operators are made here since the observations will change each new dataset (and
%           we want to run multiply minimizations for bootstrapping). 
%       data -> subsampled data in (RxSx)*(RySy) domain
%       vec  -> function that vectorizes input

%Output: 
%       out -> Data structure that contains the low-rank minimization tensor with interpolation, 
%               time to convergence, and objective value. Note that this may be
%               residuals or travel times. Also outputs the time that the routine took to
%               run

tic; 
%initialize algorithm parameters
suppress_out = params.sup_out; 
k = params.k; 
eta = params.eta; 
ns = params.ns; 
Rx = params.Rx; 
Ry = params.Ry; 
interp = params.interp;
norm_type = params.norm_type;
sigma = params.sigma; 
lambda = params.lambda; 
converged = params.converged; 
iter_crit = params.iter_crit; 
etafact = params.etafact; 
printevery = params.printevery; 
stop_crit = params.stop_crit;
g = params.g; 
I = params.I; 



%initialize data parameters
obs = data.obs; 
no_obs = data.no_obs; 
sub_tt = data.sub.tt; 
sub_res = data.sub.res; 
sub_u = data.sub.u;

%initialize
w = randn(Rx*sqrt(ns),Ry*sqrt(ns));
err          = norm(w,2);
[U,S,V]      = rsvd(w,k);
L            = U*sqrt(S);
R            = V*sqrt(S);
An = opRestriction(ns*Rx*Ry,no_obs);
Ao = opRestriction(Rx*Ry*ns,obs);

%Include uncertainties
% U = spdiags(1./sub_u(:),0, Rx*Ry*ns, Rx*Ry*ns); 
Ui = sub_u(:);
Ui(Ui==0) = 1; 
U = spdiags(1./Ui, 0, Rx*Ry*ns, Rx*Ry*ns); 
%note that here we want to keep everything else in (RxRy)x(SxSy)
%space until the end since that's what the laplacian is created is

switch(interp)
   %pick what you are matching; also scale by reported uncertainties
   case{'tt'}
        b  = Ao*(sub_tt(:));
   case{'res'}
        b  = Ao*(sub_res(:)); 
end


if strcmp(params.norm_type, 'cb')||strcmp(params.norm_type, 'lap')
    LaptLap = params.LaptLap;
    Ao_w = Ao*I; %pre-form A'A (which is weighted)
    AtA = Ao_w'*Ao_w; 
    Atb = Ao'*b;
    Lbig = params.Lbig;
    if sigma>0
%         B = Ao_w*U;
        B = Ao_w;
    end
end

obj_value = sum(svd(L*R')); 
%% run the algorithm
for i = 1:iter_crit
    count = 0;
    %prefactor everything you need (since occasionally you increase eta
    if sigma>0 &&  strcmp(norm_type, 'cb')
        A_root = [sqrt(g/2)*Lbig; sqrt(eta/2)*I]; 
%         AtA_root = A_root'*A_root;
%         R_root = chol(AtA_root); 
    end
    while err>converged && count<stop_crit
        %update observed and non-observed of w
        if strcmp(norm_type, 'lap')
            %only applies laplacian, updates only w
                wold = w(:);
%                 Q = sqrt(g/2)*Lbig; 
%                 q = -sqrt(g/2)*Lbig*Atb; 
                Q = Lbig; 
                q =-Lbig*Atb;
                [lambda,zr] = root_finder(Q, q,sigma,Ao_w, 0); 
                w = (zr+Atb); 
%                 q = [zeros(size(Ao,2),1); b]; 
%                 Q = [g*LaptLap/50000, Ao'; Ao, sparse(zeros(size(Ao,1)))]; 
%                 w = Q\q; 
%                 w = w(1:size(LaptLap,1));
%                 temp = chol(LaptLap + lambda*AtA); 
%                 w = (temp\(temp'\(-sqrt(lambda)*Atb)));
%                 lambda = 1.5*lambda; 
                err = norm(w - wold(:)); 
                werr = 0; 
                
        else
            Lold      = L;
            temp = chol(eye(k) + eta*(R'*R));
            L         = (temp\(temp'\(eta*R'*w')))';
            Rold      = R;
            temp = chol(eye(k) + eta*(L'*L)); 
            R         = (temp\(temp'\(eta*L'*w)))'; 
            LRt = vec(L*R');
            w = vec(w); 
            wold      = w;
            werr = norm(Ao*w-Ao*LRt);
            err = norm(Rold(:) - R(:)) + norm(Lold(:)-L(:)) + norm(wold(:) - w(:));
        end
        switch(norm_type)
            case{'cb'}
                %build the system: 
                if sigma==0
                    CA = eta*I + g*(LaptLap) + lambda*(AtA);
                    CA = chol(CA);
                    w = CA\(CA'\(eta*LRt+lambda*Atb)); 
                    lambda = 1.05*lambda; 
                else
                    q = -[sqrt(g/2)*Lbig*Atb; sqrt(eta/2)*(Atb - LRt)]; 
                    [lambda, z] = root_finder(A_root, q, sigma, B, lambda); 
                    w = z + Atb; 
                    
                end
            case{'l2'}
                w(no_obs) = An*LRt;
                v = Ao*LRt;
                if(norm(v-b) <= sigma)
                    w(obs) = v;
                else
                    w(obs) = sigma*(v-b)/norm(v-b) + b;
                end
                
            case{'l1'} 
                w(no_obs) = An*LRt;
                v = Ao*LRt;
                if(norm(v-b,1) <= sigma)
                    w(obs) = v;
                else
                    w(obs) = oneProjector(v-b,1,sigma) + b;
                end
            case{'lap'}
                
            otherwise
                error('unknown region type');
        end
        feas = norm(Ao*w - b) - sigma;
        if(mod(count, printevery)==1) && ~suppress_out
            fprintf('iter: %d, w-lr: %7.3f, feas: %7.3e, err: %7.3e\n', count, werr, feas, err);
        end
        %err = [err,norm(w,2)];
        obj_value = [obj_value, norm(L*R' - data.model.res)];% + norm(Lbig*vec(L*R'),2)^2]; 
        w = reshape(w, Rx*sqrt(ns),Ry*sqrt(ns)); 
        count = count+1;
%         etafact = norm(w)/feas; 
%         etafact = feas/err; (works but only until i=10 or 11
    end 
    eta = eta*etafact; 

end
out.time = toc; 
switch(params.norm_type)
    case{'lap', 'cb'}
        out.XM = w; %again ignore this for now
    case{ 'l2', 'l1'} %use LR'
        out.XM = reshape(L*R', Rx*sqrt(ns), Ry*sqrt(ns));    
end
 out.obj_val = obj_value; 
 
end
