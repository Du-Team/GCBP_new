function [link_core_nodes_indices unactivate_points] = link_core(dist_node, node_ind_point_cell, activate_node_sparse, new_unactivate_nodes,density)
% Determine the core of connections between a set of inactive nodes and activated nodes and find the inactivation point

    no_node_sparse = size(dist_node,1);
    
    unactivate_points = [];
    dist = dist_node(:,new_unactivate_nodes).*activate_node_sparse;
    dist(dist==0) = Inf;
    minDist = min(dist);
    
    link_core_nodes_indices = zeros(size(new_unactivate_nodes));
    for i = 1:length(new_unactivate_nodes)
        temp = node_ind_point_cell{new_unactivate_nodes(i)};
        unactivate_points = union(unactivate_points,temp);  
        
        minDistIndx = find(dist(:,i) == minDist(i));
        [~, secInd] = max(density(minDistIndx));
        link_core_nodes_indices(i) = minDistIndx(secInd);
    end
end