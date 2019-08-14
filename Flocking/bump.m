function out = bump(z)
    if (z >= 0) && (z < .2)
        out = 1;
    elseif (z >= .2) && (z <= 1)
        out = .5*(1 + cos(pi*((z - .2) / .8)));
    else
        out = 0;
    end
end