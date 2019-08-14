function plot_convergance(history)
    for k = 1:30
        for i = 1:25
            for j = 1:25
                figure(1)
                plot(history{k}(:,i,j))
                hold on
            end
        end
    end
    
    for k = 1:size(history{1},2)
        for i = 1:grid_size
            for j = 1:grid_size
                x = history{:,1}(1,1)(1)
            end
        end
    end

    %set initial measurement matrix
nodes_va = zeros(num_nodes,grid_size,grid_size);
observance = diff(:) < r_s(:);
for i = 1:num_nodes    
    if observance(i) %acts as observance matrix
        nodes_va(i) = F + normrnd(0,Var(i)); %sets noisy measurement
    else
        nodes_va(i) = zeros(grid_size,grid_size,1); %sets to zero
    end
end
nodes_va0 = nodes_va; %save the initial measurement