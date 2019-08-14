function plot_graph(is_2D, num_nodes, nodes, Nei_agent)
    if is_2D == true
        hold on
        subplot(2,3,1), %for 2D
        plot(nodes(:,1), nodes(:,2), 'o', 'LineWidth' ,.2 , 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'm', 'MarkerSize', 5)
        hold off
        
        for i = 1:num_nodes %graph lines between the nodes
            tmp=nodes(Nei_agent{i},:);
            for j = 1:size(nodes(Nei_agent{i},1))
                line([nodes(i,1),tmp(j,1)],[nodes(i,2),tmp(j,2)])
            end
        end
    else
        hold on
        subplot(2,3,2), %for 3D
        plot(nodes(:,1), nodes(:,2), 'o', 'LineWidth' ,.2 , 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'm', 'MarkerSize', 5)
        hold off
        
        for i = 1:num_nodes %graph lines between the nodes
            tmp=nodes(Nei_agent{i},:);
            for j = 1:size(nodes(Nei_agent{i},1))
                line([nodes(i,1),tmp(j,1)],[nodes(i,2),tmp(j,2)])
            end
        end
    end
end