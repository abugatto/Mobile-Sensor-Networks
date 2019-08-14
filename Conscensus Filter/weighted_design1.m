function weight = weighted_design1(i ,j ,Var ,Nei_agent ,Cv ,r ,r_s)
    C1w = generate_C1w(r ,r_s ,i ,j ,size(Nei_agent{i},1) ,Cv);
    
    %if i~=j and j is a member if neighbor set N_i
    if i ~= j && ismember(j, Nei_agent{i})
        weight = C1w / (Var(i) + Var(j));
    elseif i == j %for self weight Wii
        weight_sum = 0;
        for k = 1:size(Nei_agent{i},1)
            weight_sum = weight_sum + weighted_design1(i ,Nei_agent{i}(k) ,Var ,Nei_agent ,Cv ,r ,r_s);
        end
        
        weight = 1 - weight_sum;
    else %for not in neighbor set or self weight 
        weight = 0;
    end
end