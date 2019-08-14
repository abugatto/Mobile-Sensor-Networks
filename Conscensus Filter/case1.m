%===========Algorithm Parameters==============

               ALG_NUM = 2;
               
               TEST = true;
               
%===========Set Parameters====================
n = 2; %number of dimensions
F = 50; %actual value 
Cv = .01;

if ALG_NUM == 5 || ALG_NUM == 6
    num_nodes = 50; %Randomly generate nodes
    grid_size = 20; %20x20
    r = 17; %set communication range 
    r_s = 17*ones(num_nodes,1); %measurement radius
else
    num_nodes = 20; %Randomly generate nodes
    grid_size = 4; %4x4
    r = 2; %set communication range 
    r_s = 3*ones(num_nodes,1); %measurement radius
end

%Place F Value on Graph
subplot(1,3,1);
plot(.5*grid_size, .5*grid_size, 's', 'LineWidth' ,.1 , 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r', 'MarkerSize', 10)
hold on

%compute graph 
nodes = grid_size.*rand(num_nodes, n); %Randomly generate initial positions of MSN 
[Nei_agent, A] = findneighbors(nodes, n, r);
subplot(1,3,1), plot(nodes(:,1), nodes(:,2), 'o', 'LineWidth' ,.2 , 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'm', 'MarkerSize', 5)
hold on

for node_i = 1:num_nodes %graph lines between the nodes
    tmp=nodes(Nei_agent{node_i},:);
	for j = 1:size(nodes(Nei_agent{node_i},1))
        line([nodes(node_i,1),tmp(j,1)],[nodes(node_i,2),tmp(j,2)])
    end
end

%compute grid average 
x_ave = 0; 
y_ave = 0;
for i = 1:num_nodes
    x_ave = x_ave + nodes(i,1);
    y_ave = y_ave + nodes(i,2);
end
q_ave = (1/num_nodes)*[x_ave y_ave];

%Compute Variance Matrix (works for all cases in static system
Var = zeros(num_nodes,1);
for i = 1:num_nodes
    diff = norm(nodes(i,:) - q_ave);
    Var(i) = (diff^2 + Cv) / (r_s(i)^2);
end

%initial value
nodes_va = F.*ones(num_nodes,1) + normrnd(0,Var(:)); %Add measurement for each node: yi= theta + v_i
nodes_va0 = nodes_va; %save the initial measurement

%Define History
history = []; %add nodes_va_next after each iteration 
history(1,:) = nodes_va0; %each column is a history vector

%Implement Algorithm
iteration = 2;
nodes_va_next = zeros(num_nodes,1);
nodes_va_f = zeros(num_nodes,1);
while(1) %loop until convergance (then break)
    for i = 1:num_nodes %sum the neighbor weights
        %compute Wii
        if ALG_NUM == 1 %weighted alg1
            Wii = weighted_design1(i ,i ,Var ,Nei_agent ,Cv ,r ,r_s);
        elseif ALG_NUM == 2 %weighted alg2
            Wii = weighted_design2(i ,i ,Var ,Nei_agent ,Cv ,r ,r_s);
        elseif ALG_NUM == 3 || ALG_NUM == 5 %metropolis
            Wii = weighted_metropolis(i, i, Nei_agent);
        elseif ALG_NUM == 4 || ALG_NUM == 6 %max degree
            Wii = weighted_max_degree(i, i, Nei_agent);        
        end
        
        %compute sum of Wij*X_j
        sum = 0;
        for j = 1:size(Nei_agent{i},1) %iterates through neighbors
            if ALG_NUM == 1 %weighted alg1 -> (Nei_agent{i}(j) is index of j node)
                sum = sum + (nodes_va(Nei_agent{i}(j)) * weighted_design1(i ,Nei_agent{i}(j) ,Var ,Nei_agent ,Cv ,r ,r_s));
            elseif ALG_NUM == 2 %weighted alg2
                sum = sum + (nodes_va(Nei_agent{i}(j)) * weighted_design2(i ,Nei_agent{i}(j) ,Var ,Nei_agent ,Cv ,r ,r_s));
            elseif ALG_NUM == 3 || ALG_NUM == 5 %metropolis
                sum = sum + (nodes_va(Nei_agent{i}(j)) * weighted_metropolis(i, Nei_agent{i}(j), Nei_agent));
            elseif ALG_NUM == 4 || ALG_NUM == 6 %max degree
                sum = sum + (nodes_va(Nei_agent{i}(j)) * weighted_max_degree(i, Nei_agent{i}(j), Nei_agent));
            end
        end
        
        %get estimates X_i
        nodes_va_next(i) = Wii*nodes_va(i) + sum;
    end
    
    convergance = false; 
    conv_val = nodes_va_next(1); 
    for i = 2:num_nodes %check if converged within .001%
        if (conv_val - nodes_va_next(i)) < .001 %breaks if not in conscensus
            if i == 10
                convergance = true;
                nodes_va_f = conv_val; %set final converged estimate
            end
        end
    end
    
    if convergance == true
        break
    end
    
    history(iteration,:) = nodes_va_next; %add estimate to history
    nodes_va = nodes_va_next; %iterate estimates
    iteration = iteration + 1;
    
    if TEST == true
        plot_data(TEST, num_nodes, history, nodes_va_next, nodes_va_f, nodes_va0)
    end
end

if TEST == false
    plot_data(TEST, num_nodes, history, nodes_va_next, nodes_va_f, nodes_va0);
end