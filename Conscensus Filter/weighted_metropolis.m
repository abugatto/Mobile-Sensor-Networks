function weight = weighted_metropolis(i, j, Nei_agent)
    if i ~= j && ismember(j, Nei_agent{i})
        weight = 1 / (1 + max(size(Nei_agent{i},1), size(Nei_agent{i},1)));
    elseif i == j
        weight_sum = 0;
        for k = 1:size(Nei_agent{i}, 1)
            weight_sum = weight_sum + weighted_metropolis(i, Nei_agent{i}(k), Nei_agent);
        end
        
        weight = 1 - weight_sum;
    else
        weight = 0;
    end 
end