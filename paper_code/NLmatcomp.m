function [f, g, o] = NLmatcomp(x,Ao,params)
Rx = params.Rx; 
Ry = params.Ry; 
ns = params.ns; 
b = params.b;
k = params.k; 
nab = params.Lbig; 
lambda = params.lambda;
gamma = params.gamma; 
AtA = params.AtA;
LtL = params.LtL; 
e = Rx*sqrt(ns)*k;%unpack all entries up to r
l = x(1:e);%proxy for l gradient
r = x(e+1:end);%proxy for r gradient
L = reshape(l,Rx*sqrt(ns),k); %is this right? 
R = reshape(r,k, Ry*sqrt(ns)); %run by sasha
LR = L*R; 
o = norm(LR - params.true); 
% if isempty(g)
%     f1 = params.afun(L*R');
%     f2 = 0;
% else 
%     fp = params.afunT(g);
%     f1 = [vec(fp*R); vec(fp'*L)];
%     f2 = vec(fp);
% end
%function evaluation
f = .5*norm(L, 'fro')^2 + .5*norm(R, 'fro')^2+...
    gamma*norm(nab*LR(:),2)^2+...
    lambda*norm(Ao*LR(:) - b,2)^2; 
% solve for gradient
% C = nab*(kron(R', speye(size(L,1)))); 
% D = Ao*kron(R', speye(size(L,1))); 
%l gradient
% l = l + 2*gamma*C'*(C*l)+2*lambda*D'*(D*l) - 2*lambda*D'*b; 
temp = 2*gamma*reshape(LtL*LR(:), Rx*sqrt(ns), Ry*sqrt(ns))*R';
temp2 = 2*lambda*reshape(AtA*LR(:), Rx*sqrt(ns), Ry*sqrt(ns))*R';
temp3 =  2*lambda*reshape(Ao'*b, Rx*sqrt(ns), Ry*sqrt(ns))*R'; 
l = l + temp(:) + temp2(:) - temp3(:); 

%repeat for R
% C = nab*(kron(speye(size(R,2)),L)); 
% D = Ao*kron(speye(size(R,2)),L); 
% r = r + 2*gamma*C'*(C*r) + 2*lambda*D'*(D*r)- 2*lambda*D'*b;
C = kron(speye(size(R,2)), L'); 
r = r + 2*gamma*C*(LtL*LR(:)) + 2*lambda*C*(AtA*LR(:)) - 2*lambda*C*(Ao'*b);
g = [l; r]; 

end