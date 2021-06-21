function test_suite = test_oplWindow1D
%test_oplWindow1D  Unit tests for 1d windowing operators
initTestSuite;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
function test_oplWindow1D_fd
%%
   for n=30:3:42
       for p=2:3
           for h=0:3
	       A=oplWindow1Dfd(n,p,h);
               utest(A);
	   end
       end
   end
end

function test_oplWindow1D_tpr
%%
   for n=30:3:42
       for p=2:3
           for h=0:3
	       A=oplWindow1Dtpr(n,p,h);
               utest(A);
	   end
       end
   end
end

function test_oplWindow1D_avg
%%
   for n=30:3:42
       for p=2:3
           for h=0:3
	       A=oplWindow1Davg(n,p,h);
               utest(A);
	   end
       end
   end
end

