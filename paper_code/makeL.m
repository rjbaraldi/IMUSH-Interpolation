function [Lbig,Lrhs] = makeL( Rx,Ry,evind )
[~, ~, L]      = make_laplacian(Rx,5,Ry,5);%is this right?
Lbig           = sparse([],[],[],length(evind)*size(L,1),length(evind)*size(L,1),0); % used in inversion, Laplacian operators to smooth model for each event
Lrhs           = [];
Lr             = 1;
Lc             = 1;

for k = 1:length(evind) % loop over events
    Lbig([1:size(L,1)]+Lr-1,[1:size(L,2)]+Lc-1)=L;
    Lr   = Lr+size(L,1);
    Lc   = Lc+size(L,2);
    Lrhs =[Lrhs;zeros(size(L,2),1)];  %%%% add to rhs
end

end

