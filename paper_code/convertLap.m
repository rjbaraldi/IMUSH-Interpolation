function newL = convertLap(L, Rx,Ry, ns, sw_cond)
%not quite sure how to use the opPermutation thing so i'll do it on my own
%for now
vec = @(x) x(:);
[m,n] = size(L); 
newL = spalloc(m,n, numel(find(L~=0)));
switch(sw_cond)
    case{'full'}
        for i = 1:m
            newL(i,:) = sparse(vec(reshape(permute(reshape(full(L(i,:)), Rx,Ry,sqrt(ns),sqrt(ns)),[1 3 2 4]),Rx*sqrt(ns),Ry*sqrt(ns)))');    
        end
 case{'index'}
            newL = vec(reshape(permute(reshape(full(L), Rx,Ry,sqrt(ns),sqrt(ns)),[1 3 2 4]),Rx*sqrt(ns),Ry*sqrt(ns)))';
end   


end