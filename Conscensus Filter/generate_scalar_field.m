function F = generate_scalar_field(dim, scale, theta, corr, var, mu)
    meshdim = scale(2):scale(1):scale(3); %mesh dimensions
    [X, Y] = meshgrid(meshdim(1,:),meshdim(1,:)); %create mesh axes values
    
    phi = zeros(size(mu,1), size(X,1)*size(Y,1),1); %4x625 matrix
    for i = 1:size(mu,1)
        covxy = corr(i)*sqrt(var(i,1)*var(i,2));
        cov = [var(i,1) covxy; covxy var(i,2)];
        phi(i,:) = mvnpdf([X(:) Y(:)], mu(i,:), cov); %pdf is 1x625 vector
    end
    
    F = theta * phi; %weighted distribution mixture F = [1x4]*[4x625]
    
    %plot mesh and color
    num_cells = size(F,2);
    F = reshape(F, size(X,1), size(Y,1)); %form 25x25 grid using gaussian mixture
    
    subplot(2,3,1);
    pcolor(meshdim, meshdim, -F)
    hold on
    
    subplot(2,3,2);
    surf(meshdim, meshdim, -F)
    hold on
end

