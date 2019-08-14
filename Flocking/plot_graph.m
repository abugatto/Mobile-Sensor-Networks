function plot_graph(N,r,pos,nlist)
    figure
    hold on
    
    %use imaginary numbers as tool to draw line segments
    for k = 1:max(size(nlist)) %iterates through neighbor list
        m = nlist(1,k); %all indices (node 1,2,3,...etc)
        n = nlist(2,k); %all node neighbors (node1.1 1.2 1.3 ...etc)
        %plot all line segments
        plot([pos(m,1) + 1i*pos(m,2), pos(n,1) + 1i*pos(n,2)], 'b', 'LineWidth', 1);
    end
    
    %plot all points onto the complex space
    plot(pos(:,1) + 1i*pos(:,2), '>', 'Color', 'r', 'LineWidth', 2, 'MarkerSize', 5);
    hold on;
    
    %{
    for w = 1:N
        theta = 0:pi/50:2*pi; %plot communication area for nodes
        x = r*cos(theta) + pos(w,1);
        y = r*sin(theta) + pos(w,2);
        plot(x, y, 'Color', 'g', 'LineWidth', 1);
        hold on;
    end
    hold off;
    %}
    
    h = get(gcf, 'CurrentAxes'); %get current figure axes
    set(h, 'FontName', 'Helvetica', 'FontSize', 10);
end


