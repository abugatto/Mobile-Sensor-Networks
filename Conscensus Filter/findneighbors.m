function [Nei_agent, A] = findneighbors(nodes, n, r)

%This function is to find alpha and beta neighbors
%Created by Anthony Bugatto

% Inputs:  positions of nodes (nodes),
          %dinensions of the simulation(n)
          %active range for alpha agents(r)
          
% Outputs: indices of alpha neighbors (Nei_agent)
           %Adjacency Matrix (A)

%*******Find neighbors of alpha agent*************

num_nodes = size(nodes,1);
dif = cell(num_nodes,1); % save the difference between each alpha agent and all other nodes 
                         % each element of cell is a matrix(size:num_nodes x n)
distance_alpha = zeros(num_nodes,num_nodes); % save the distance (norm) between each agent and all other nodes
                                % each column for one node      
Nei_agent = cell(num_nodes,1); %Save the indices of neighbors of each agent 
                               %each element of cell is a matrix (maximum size num_nodes x 1)
                                
for i = 1:num_nodes
    dif{i} = repmat(nodes(i,:),num_nodes,1) - nodes;
    tmp = dif{i}; %recall cell i th of dif
    
    for j = 1:num_nodes
        d_tmp(j,:) = norm(tmp(j,:)); %compute distance between each alpha agent and all other nodes
    end
    
    distance_alpha(i,:)= d_tmp; 
end

for k = 1:num_nodes
    Nei_agent{k} = find(distance_alpha(:,k) < r & distance_alpha(:,k) ~= 0); %find the neighbors of agent i
end

%+++++++++++++++++BUILDING THE ADJACENCY MATRIX+++++++++++++++++++++

A = zeros(num_nodes,num_nodes);

for i = 1:num_nodes
    for j = 1:num_nodes
        if i ~= j
            dist_2nodes = norm(nodes(j,:) - nodes(i,:));
            if dist_2nodes < r && dist_2nodes ~= 0
                A(i,j) = 1; 
            end
        end
    end 
end


