function [ normBx,x, df] = mynorm( lam, AtA, B, Atq, sig2 )

% R = chol(AtA + lam*(B'*B));
% x  = R\(R'\Atq);
% dx = 1; 
%x = (AtA+lam*(B'*B))\Atq; 
[x, ~] = pcg(AtA+lam*(B'*B), Atq, 1e-8, 100);
% dx = -B'*(B)*(AtA + lam*(B'*B))^(-2)*Atq; % \nabla x(lam)
normBx = norm(B*x)-sig2;
dx = imag(x)/imag(lam); %get gradient of x(lam)
df = real((B'*(B*x))'*dx/norm(B*x)); %chain rule to get gradient of f(lam)
%  fprintf('    Rf iter:   Lambda: %d    Compare to 1e-6: %d \n', lam, abs(normBx));
end

