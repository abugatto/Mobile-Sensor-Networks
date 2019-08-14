function signorm = sigma_norm(in)
    epsilon = .1; %needs to be less than 1
    signorm = (1/epsilon)*(sqrt(1 + epsilon*(norm(in,2)^2)) -1);
end

%test and comparison
%s1 = sigma_norm([3,5],[7,5]) = 6.1245
%norm = 4

%s2 = sigma_norm([2,8],[4,5]) = 5.1658
%norm = 3.6056


