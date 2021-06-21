function test_suite = test_opLBFGS
%test_opLBFGS  Unit tests for opLBFGS.
initTestSuite;
end

function test_opLBFGS_symmetric

   rng('default');

   n = 10;
   mem = 5;
   npairs = 7;

   % Set up forward and inverse LBFGS operators
   B = opLBFGS(n, mem);
   B.update_forward;
   for ~ = 1 : npairs
     s = rand(n, 1);
     y = rand(n, 1);
     B = update(B, s, y);
   end

   % Check that B and H are symmetric
   % (cheap test taken from MINRES)
   y  = rand(n, 1);
   w  = B \ y;
   r2 = B \ w;
   s  = w' * w;
   t  = y' * r2;
   z  = abs(s - t);
   epsa = (s + eps) * eps^(1/3);
   assert(z <= epsa);

   w = B * y;
   r2 = B * w;
   s  = w' * w;
   t  = y' * r2;
   z  = abs(s - t);
   epsa = (s + eps) * eps^(1/3);
   assert(z <= epsa);

end

function test_opLBFGS_definite

   rng('default');

   n = 10;
   mem = 5;
   npairs = 7;

   % Set up forward and inverse LBFGS operators
   B = opLBFGS(n, mem);
   B.update_forward
   for ~ = 1 : npairs
     s = rand(n, 1);
     y = rand(n, 1);
     B = update(B, s, y);
   end

   % Check that B and H are positive definite
   % (cheap test taken from MINRES)
   v = rand(n, 1);
   Bv = B \ v;
   vBv = v' * Bv;
   assert(vBv > 0);

   Bv = B * v;
   vBc = v' * Bv;
   assert(vBv > 0);

end

function test_opLBFGS_inverse

   rng('default');

   n = 10;
   mem = 5;
   npairs = 7;

   % Set up forward and inverse LBFGS operators
   B = opLBFGS(n, mem);
   B.update_forward
   for ~ = 1 : npairs
     s = rand(n, 1);
     y = rand(n, 1);
     B = update(B, s, y);
   end

   Bmat = B * eye(n);
   assertElementsAlmostEqual(B \ Bmat, eye(n));
   Bmat = B \ eye(n);
   assertElementsAlmostEqual(B * Bmat, eye(n));

end
