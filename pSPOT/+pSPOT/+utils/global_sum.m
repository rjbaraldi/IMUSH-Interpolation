function y = global_sum(y)
% loop through labs checking if its ready to send data, and receive if it
% is

if labindex == 1   %sum all results on lab 1
    labs = 2:numlabs;
    while ~isempty(labs)
        if( labProbe(labs(1)) )
            y = y + labReceive(labs(1));
            labs(1) = [];
        else
            labs = circshift( labs, [0 -1] );
        end
    end
else % Other labs
    labSend(y,1);
end