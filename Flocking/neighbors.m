function [edges, nlist] = neighbors(N,pos,r)
    T = (N*(N-1))/2; %total possible edges
    
    num_edges = 0;
    nlist = zeros(2,T);
    
    for p = 1:N
        for q = p+1:N
            if(norm(pos(p,:) - pos(q,:), 2) <= r) %computing euclidean norm
                num_edges = num_edges + 1; 
                
                %set nlist value equal to the pair of neighbor vertices
                nlist(:,num_edges) = [p,q]'; 
            end
        end
    end
    
    nlist = nlist(:,1:num_edges); %cutting off all the left over terms
    edges = sparse(nlist);
    
    plot_graph(N, r, pos, nlist);
end