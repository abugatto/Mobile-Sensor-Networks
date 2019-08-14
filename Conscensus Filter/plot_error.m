function  plot_error(F, nodes_va_next, nodes_va_f, scale)
    meshdim = scale(2):scale(1):scale(3); %mesh dimensions
    
    error = zeros(size(meshdim,2), size(meshdim,2));
    final = zeros(size(meshdim,2), size(meshdim,2));
    for k = 1:size(meshdim,2)
        for l = 1:size(meshdim,2)
            %final(k,l) = nodes_va0(1,k,l);
            final(k,l) = nodes_va_next(10,k,l); %assume one node is correct
            error(k,l) = F(k,l) - final(k,l);
        end
    end
    
    hold on
    subplot(2,3,4) %new grid
    pcolor(meshdim, meshdim, -final)
    hold off

    hold on
    subplot(2,3,5) %new mesh
    surf(meshdim, meshdim, -final)
    hold off

    hold on
	subplot(2,3,3) %stem error
    stem(-error(:))
    hold off
    
    hold on
    subplot(2,3,6) %mesh error
    surf(meshdim, meshdim, -error)
    hold off
end

%{
subplot(1,3,1), surf(meshdim, meshdim, -final)
subplot(1,3,2), surf(meshdim, meshdim, -F)
subplot(1,3,3), surf(meshdim, meshdim, F-final)
%}
