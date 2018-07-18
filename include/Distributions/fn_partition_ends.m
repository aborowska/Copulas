function [s1, s2] = fn_partition_ends(partition, d, s)
    s1 = partition(s);
    S = length(partition);
    if s < S
        s2 = partition(s+1)-1;
    else
        s2 = d;
    end
end