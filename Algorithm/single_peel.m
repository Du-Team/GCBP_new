function [activate_node_sparse old_unactivate_nodes new_unactivate_nodes] =single_peel(iters, precentile, node_sparse, density, activate_node_sparse, old_unactivate_nodes);
% Strip data points

% Density ordering
    dens_mat = [node_sparse,density];
    [~,orderByDensity] = sort(density);
    nodesSortedByDens = dens_mat(orderByDensity,:);
    
    no_node_sparse = size(nodesSortedByDens,1);
    
    % density threshold
    index_prcentile = ceil(no_node_sparse* (1-(1-precentile)^iters)); % Find the index value corresponding to 10%
    threshold_value = nodesSortedByDens(index_prcentile,end); % Find the b-value corresponding to this index as the threshold tau
    
    % Node activation
    unactivate_nodes = find(density <= threshold_value)';
    new_unactivate_nodes = setdiff(unactivate_nodes, old_unactivate_nodes);
    old_unactivate_nodes = unactivate_nodes;
    activate_node_sparse(new_unactivate_nodes)=0;
end
