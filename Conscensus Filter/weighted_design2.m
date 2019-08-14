function weight= weighted_design2(i ,j ,Var ,Nei_agent ,Cv ,r ,r_s)
    C2w = (.5*Cv) / (r_s(i)^2);
    
    %if i~=j and j is a member if neighbor set N_i
    if i ~= j && ismember(j, Nei_agent{i})
        weight = (1 - weighted_design2(i ,i ,Var ,Nei_agent ,Cv ,r ,r_s)) / size(Nei_agent{i},1);
    elseif i == j %for self weight Wii
        weight = C2w / Var(i);
    else %for not in neighbor set or self weight 
        weight = 0;
    end
end 