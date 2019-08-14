function C1w = generate_C1w(r ,r_s ,i ,j ,num_nei ,Cv) %weight parameter
    %if communication and sensing distances are equal
	if (r_s(i) == r_s(j)) && (r_s(i) == r) && (r_s(j) == r)
        C1w = 2*(Cv / ((r_s(i)^2) * length(num_nei-1))); %set legnth as N-1 for inequality
    else
        C1w = ((1/length(num_nei-1))) * (Cv/(r_s(i)^2) + Cv/(r_s(j)^2));
	end   
end