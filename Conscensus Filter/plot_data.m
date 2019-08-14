function plot_data(TEST, num_nodes, history, nodes_va_next, nodes_va_f, nodes_va0)
    %Plot Convergance
    subplot(1,3,2);
    for i = 1:num_nodes
        plot(history(:,i) ,'color' ,rand(1,3))
        hold on
    end

    %Plot Final Estimates vs Initial Estimates
    if TEST == true
        final = nodes_va_next; %TEST CODE: analize before convergance
    else
        final = nodes_va_f*ones(10,1); %create array of conscensus values for plotting
    end

    subplot(1,3,3);
    plot(final,'ro')
    hold on
    plot(final,'r')
    hold on
    plot(nodes_va0,'bo')
    hold on
    plot(nodes_va0,'b')
    drawnow
    hold off
end