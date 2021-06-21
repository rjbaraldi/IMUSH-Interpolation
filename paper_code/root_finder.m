function [ lam, z ] = root_finder(A, q, sig2, B, lam )
% returns z and lam in  ||Az - q||^2 + lam ||Bz||^2 
%   that ensures ||Bz|| = sig2
%   Detailed explanation goes here (fill in)
Atq = A'*q;
AtA = A'*A; 
%R = chol(AtA); 
h = 1e-8; 
tol = h; 
converged = 0; 
iter = 0; 
[m,n] = size(A);
%zunc = R\(R'\Atq);
% zunc = lsqr(m,n,A, q, 0, 1e-8, 1e-8, 100, 100,0);
zunc = lsqr(A, q, 1e-8, 100);
if(norm(B*zunc) <= sig2)
    z = zunc; 
    lam = 0.0; 
else
    myF = @(lam)mynorm( lam, AtA, B, Atq, sig2);
%    [lam,~] = fzero(myF, lam);
    while(~converged && iter<100)
        iter = iter + 1;
        [fv,z] = myF(lam);
        [~, ~, fd] = myF(lam+1i*h); 
        lam = max(lam - fv/fd, 0);
        % if lam = 0, something is wrong. print out a warning and exit 
        if lam==0
            disp('lam=0: Root-finding problem is infeasible');
            break
        end
        converged = norm(B*z)-sig2 < tol;
        fprintf('iter: %d, lam: %7.3e, err: %7.3e\n', iter, lam, norm(B*z)-sig2);
    end

    [~, z, ~] = myF(lam);
end
%     iter = 0; 

%(chain rule of lambda)
% as B'(Bx(lam))/||Bx(lam)|| \nabla x(lam)
% try fzero
% and if you get a reliable result with fzero
% but not yours
% then use it to get gradient of ||Bx(lam)||_2
% then use that to implement rootfinder for ||Bx(\lam)||_2 = sig2
% 
% converged = 0;
% % lam = 13; 
% tol = 1e-6;
% iter = 0;
% lam = 0;
% while(~converged) %changed (usually kept to less than 3)
%     iter = iter + 1;
%     R = chol(AtA + lam*(B'*B)); %the condition # here blows up after 1 iteration
% %     R = chol(AtA + lam*(speye(size(A, 2)))); 
%     z = R\(R'\Atq); 
%     ql = R'\z;
% %     lam = max(lam + (norm(z)/norm(ql))^2 *(norm(z)-sig2)/sig2, 0);  
% %     converged = abs(norm(z)-sig2) < tol;
% 
% 
% % SMART WAY:
%     lam = max(lam + (norm(B*z)/norm(ql))^2 *(norm(B*z)-sig2)/sig2, 0);  
%     converged = abs(norm(B*z)-sig2) < tol; 
% %     if ~sup_out
%         fprintf('    Rf iter #: %d   Lambda: %d    Compare to 1e-6: %d \n', iter, lam, abs(norm(B*z)-sig2));
% %     end
% end
% fprintf('----------- Root finding complete for this w --------\n');
%iter

    