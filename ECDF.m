function p = ECDF(y)
    [~, ind] = sort(y);
    [~,p] = sort(ind);
    p = p/(T+1);
end