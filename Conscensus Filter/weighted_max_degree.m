function weight = weighted_max_degree(i, j, Nei_agent)
    if i ~= j && ismember(j, Nei_agent{i})
        weight = 1 / size(Nei_agent,1);
    elseif i == j
        weight = 1 - (size(Nei_agent{i},1) / size(Nei_agent,1));
    else
        weight = 0;
    end 
end