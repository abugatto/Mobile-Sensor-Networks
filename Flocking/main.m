%=======================================================
%Anthony Bugatto
%CS 455: Mobile Sensor Netorks
%Project 1: Flocking
%=======================================================

                    clc,clear
                    close all
                    
%================ALGORITHM NUMBER=======================

                    ALG_NUM = 4;    
                    
%=============PARAMETER OF SIMULATION===================

d = 15;%Set desired distance among sensor nodes
k_scale = 1.2;%Set the scale of MSN
r = k_scale*d; %Set the active range 
r_prime = .22*k_scale*r; %Set the active range of beta agent
epsilon = 0.1; %Set a constant for sigma norm
num_nodes = 100; %Set number of sensor nodes
n=2; %Set number of dimensions

grid_size = 0;
if ALG_NUM ~= 1 && ALG_NUM ~= 2
    grid_size = 150; 
else
    grid_size = 50; 
end

nodes = grid_size.*rand(num_nodes, n) + grid_size.*repmat([0 1], num_nodes, 1); %Randomly generate initial positions of MSN 
p_nodes = zeros(num_nodes,n); %Set initial velocties of MSN
delta_t_update = 0.08; %Set time step
t = 0:delta_t_update:7;% Set simulation time

if ALG_NUM == 5 %=================SET OBSTACLES==========================
    obstacles = [50, 100; 150 80; 200, 230; 280 150]; %set positions of obstacles
    Rk = [20; 10; 15; 8]; %Radii of obstacles
    num_obstacles = size(obstacles,1); %Find number of obstacles
end

if ALG_NUM ~= 1 %=========SET A STATIC TARGET================
    qt1 = [150 150]; %Set position of the static target (gamma agent)
	pt1= [0 0]; %Set initial velocity of the target
end
    
nodes_old = nodes; %KEEP previous positions of MSN
q_mean = zeros(size(t,2), n); %Save positions of COM (Center of Mass)
p_mean = zeros(size(t,2), n); %Save velocities of COM (Center of Mass)
Connectivity = []; %save connectivity of MSN
q_nodes_all = cell(size(t,2), num_nodes); %creates cell array to store history of system pos
p_nodes_all = cell(size(t,2), num_nodes); % -   -   -   -   -   -   -   -   -   -   -   -vel

 
nFrames = 20; %set number of frames for the movie
mov(1:nFrames) = struct('cdata', [],'colormap', []); %Preallocate movie structure

for iteration = 1:length(t)
    if ALG_NUM ~= 1
        if ALG_NUM == 3
           %Sinewave Trajectory of a moving target
           qt_x1 = 50 + 50*t(iteration);
           qt_y1 = 295 - 50*sin(t(iteration));
        elseif ALG_NUM == 4
            %Circle Trajectory of a moving target
            qt_x1 = 310 - 160*cos(t(iteration)); 
            qt_y1 = 255 + 160*sin(t(iteration)); 
        elseif ALG_NUM == 5
            %Line Trajectory of a moving target
            qt_x1 = 200 + 130*t(iteration);
            qt_y1 = 200 + 1*t(iteration);
        end

        if ALG_NUM == 3 || ALG_NUM == 4 || ALG_NUM == 5 
            %compute position of target
            qt1(iteration,:) = [qt_x1, qt_y1];
        
            %compute velocities of target
            if iteration > 1
                pt1(iteration,:) = (qt1(iteration,:) - qt1(iteration-1,:)) / delta_t_update;
            else
                continue 
            end 
        end
        
        plot(qt1(:,1),qt1(:,2),'ro','LineWidth',2,'MarkerEdgeColor','r','MarkerFaceColor','r', 'MarkerSize',4.2)
        hold on
    end

    %=====================(LOOP) UPDATE PROCESS for MSN=====================
    
    if ALG_NUM == 1 %fragmentation
        [Nei_agent, A] = findneighbors1(nodes, r);
        [Ui] = inputcontrol_Algorithm1(num_nodes, nodes, Nei_agent, n, epsilon, r, d, p_nodes); 
    elseif ALG_NUM == 2 % static target
        [Nei_agent, A] = findneighbors1(nodes, r);
        [Ui] = inputcontrol_Algorithm2(ALG_NUM, num_nodes, nodes, Nei_agent, n, epsilon, r, d, qt1, pt1, p_nodes); 
    elseif ALG_NUM == 3 || ALG_NUM == 4 %goal following
        [Nei_agent, A] = findneighbors1(nodes, r);
        [Ui] = inputcontrol_Algorithm2(ALG_NUM, num_nodes, nodes, Nei_agent, n, epsilon, r, d, qt1(iteration,:), pt1(iteration,:), p_nodes); 
    elseif ALG_NUM == 5 %goal following and obstacle avoiding
        [Nei_agent, Nei_beta_agent, p_ik, q_ik, A] = findneighbors5(nodes_old,nodes, r, r_prime, obstacles, Rk, n, delta_t_update);
        [Ui] = inputcontrol_Algorithm5(num_nodes, nodes, Nei_agent, n, epsilon, r, r_prime, d, k_scale, Nei_beta_agent, p_ik, q_ik, obstacles, qt1, pt1, p_nodes);
    end
    
    p_nodes = (nodes - nodes_old)/delta_t_update; %COMPUTE velocities of sensor nodes
    p_nodes_all{iteration} = p_nodes; %SAVE VELOCITY OF ALL NODES
    nodes_old = nodes;
    nodes = nodes_old + p_nodes*delta_t_update + .5*Ui*delta_t_update*delta_t_update;
    q_mean(iteration,:) = mean(nodes); %Compute position of COM of MSN
    
    if ALG_NUM ~= 1
        plot(q_mean(:,1),q_mean(:,2),'ro','LineWidth',2,'MarkerEdgeColor','k', 'MarkerFaceColor','k','MarkerSize',4.2)
        hold on
    end
    
    p_mean(iteration,:) = mean(p_nodes); %Compute velocity of COM of MSN
    q_nodes_all{iteration} = nodes; 
    Connectivity(iteration)= (1 / num_nodes) * rank(A);
   
    if ALG_NUM == 5 %Draw obstacles
         phi = 0:.1:2*pi;
         for k = 1:num_obstacles
            X = Rk(k)*cos(phi);
            Y = Rk(k)*sin(phi); 
            plot(X+obstacles(k,1),Y+obstacles(k,2),'r',nodes(:,1),nodes(:,2), 'g>') 
            fill(X+obstacles(k,1),Y+obstacles(k,2),'r')
            axis([0 250 0 80]);
            hold on
         end
    end
     
    %================= PLOT and LINK SENSOR TOGETHER ===============
    plot(nodes(:,1),nodes(:,2), '.')
    hold on
    plot(nodes(:,1),nodes(:,2), 'k>','LineWidth',.2,'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',5) 
    hold off
    
    for node_i = 1:num_nodes
        tmp=nodes(Nei_agent{node_i},:);
        for j = 1:size(nodes(Nei_agent{node_i},1))
            line([nodes(node_i,1),tmp(j,1)],[nodes(node_i,2),tmp(j,2)])
        end
    end
    
    mov(iteration) = getframe;
    hold off
end

%==========================VIDEO SIMULATION=============================
%{
v = VideoWriter('flocking.avi', 'MPEG-4'); %Make movie
open(v)
writeVideo(v,mov);
%}
%========================PLOT VELOCITY OF MSN===========================

p_each_nodes = [];
for i = 2:size(t,2) %iterates through the timesteps for the history cell matrix
    tmp7 = p_nodes_all{i};
    for j = 1:num_nodes
        if j == 1 %Plot velociy of sensor node 1; you can change this number to plot for other nodes 
            p_each_nodes(i) = norm(tmp7(j,:));
            figure(3), plot(p_each_nodes, 'b')
            hold on
        end
    end
end

figure(4), plot(Connectivity)
grid on

%========================PLOT TRAJECTORY OF SENSOR NODES=============== 

for i = 2:length(q_nodes_all)
    tmp8 = q_nodes_all{i};
    figure(5), plot(tmp8(:,1), tmp8(:,2), 'k.') 
    hold on
end

hold on
plot(nodes(:,1), nodes(:,2), 'm>', 'LineWidth', .2, 'MarkerEdgeColor', 'm', 'MarkerFaceColor', 'm', 'MarkerSize', 5)

%========================PLOT TRAJECTORY OF COM AND TARGET===============

figure(6), plot(q_mean(:,1), q_mean(:,2),'k.') 
hold on

if ALG_NUM ~= 1 || ALG_NUM ~= 2
    plot(qt1(:,1), qt1(:,2), 'r.')
end
