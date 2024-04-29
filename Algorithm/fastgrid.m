function ptClusters = fastgrid(data,no_grid,precentile,eps,max_iters)



verbose = 0; % Whether or not to print
NdChangeThres = 1; 
% Allow two iterations of node change Take 1 to indicate termination 
% if there is no transformation before or after the two iterations, 
% and 0 to indicate that no judgment is made.


[no_objs,no_dims] = size(data); % no_objs: data number, no_dims: data dimension
activate_data = ones(no_objs,1); 
% The activation status of the data points, 
% with 1 being a core point still involved in the calculation 
% and 0 being a stripped point no longer involved in the calculation

activate_data_remat =repmat(activate_data,[1 2^no_dims])'; 
% Matrix of (2^m)*n tensors with activation of the ith data point relative 
% to the (2^m) nearest neighbor nodes accessible via activate_data_remat(:,i)


%%%%%%%%%%%%%%%%%%%%%%%%%% Data preprocessing %%%%%%%%%%%%%%%%%
 %
 %

 data=((data-min(data))./(max(data)-min(data))).*(no_grid-1)+1;% data scaling

%Fine tuning to stay within the grid
data_temp = data;
data_temp(data_temp==no_grid) = no_grid-(10^(-5));

% Generate tensor and repetitive forms of data for subsequent processing
data_tensor = permute(repmat(data_temp,...  % (2^m)*m*n tensor
    [1 1 2^no_dims]), [3 2 1]);  % The ith data repeated 2^m times can be accessed via data_tensor(:,:, i)

%
%
%%%%%%%%%%%%%%%%%%% Data preprocessing END %%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%Iterate through the neighboring nodes of each data point%%%%%%%%%
%
%  Find the node coordinates of the nodes around the data

data_floor = floor(data_temp);% data_floor is the node with the smallest coordinates of the data
% Maps 2^no_dims node coordinates for subsequent use.
data_floor_tensor = repmat(data_floor,[1 1 2^no_dims]); 
% The result is an n*m*(2^m) tensor, but that approach is not convenient for traversing its 
% surrounding nodes through the data points, so the following transformation is done
data_floor_tensor = permute(data_floor_tensor,[3 2 1]); 

% bin_sequence ordered binary array, i.e., if the dimension is 2, the
% [0 0
%  0 1
%  1 0
%  1 1]
bin_sequence = dec2bin(0:(2^no_dims-1), no_dims)-48;
bin_tensor = repmat(bin_sequence,[ 1  1 no_objs]); % (2^m)*m*n µÄtensor£»
clear bin_sequence;
% data_nodes(:,:, i) denotes the coordinates of all the nearest neighbor nodes of the ith data for each data point
data_nodes_tensor = data_floor_tensor+bin_tensor;
clear data_floor_tensor;
clear bin_tensor;
data_nodes_mat = reshape(permute(data_nodes_tensor, [1 3 2]),no_objs*2^no_dims, []);

%
%
%%%%%%%%%%%%%%%%%%Iterate through the neighboring nodes of each data point END %%%%%%%%%%

%%%%%%%%%%%%%%%%%%%  Calculate the individual gains of each node %%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
fd_tensor = 1- abs(data_tensor - data_nodes_tensor);    
single_gain_tensor = prod(fd_tensor,2);
clear fd_tensor;
single_gain_mat = reshape(single_gain_tensor,2^no_dims, no_objs);
single_gain_onedim = reshape(single_gain_mat,no_objs*2^no_dims, []); %
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Calculate the individual gains of each node END %%%%%%%%%%%%%

%%%%%%%%%%%% Calculation of density    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%

% Find the coordinates of useful nodes and set the activation array of node nodes
[node_sparse, ~, IC] = uniqueRows(data_nodes_mat);  % Real nodes for computation
no_node_sparse = size(node_sparse,1); % Number of nodes for which density exists

if length(no_node_sparse) < 1000
    min_cluster_size = 0;
else
    min_cluster_size = 0;
end

membership_remat = ind2vec(IC')';

density = sum(bsxfun(@times, membership_remat,single_gain_onedim))';
clear single_gain_onedim;
%
%
%%%%%%%%%%%%% Calculation of density  END  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%% Iterate over the data points corresponding to each node %%%%%%%%%%%%%%%%%%%%%%
% 
%

sp = ones(no_objs,1).*2^no_dims; % Spacing method for matrix chunking
membership_remat_cell = mat2cell(membership_remat,sp); % The cell where the matrix is stored after chunking
clear sp;
membership_cell = cellfun(@(x) max(x), membership_remat_cell,'UniformOutput',false); 
% Each block finds the maximum value of its columns, 
% i.e. 2^no_dims rows are merged into 1 row
membership = cell2mat(membership_cell);
sp_nodes = ones(no_node_sparse,1);
node_ind_point_cell = cellfun(@(x) find(x), mat2cell(membership',sp_nodes),'UniformOutput',false);

%
%
%%%%%%%%%%%%%%%% Iterate over the data points corresponding to each node  END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Layer-by-layer peeling operation

dist_node = squareform(pdist(node_sparse,'Chebychev'));

activate_node_sparse= ones(no_node_sparse,1);
old_unactivate_nodes = [];
cluster_uf = UF(no_node_sparse);
iters = 1;
densIters = density;
meanBorderDens = [];
prevIterNdLen = 0;

border_nodes_indices_per_iteration ={}; 
core_nodes_indices_per_iteration = {};

while(max_iters>=iters) 
    [activate_node_sparse old_unactivate_nodes new_add_unactivate_nodes] = single_peel(iters, precentile, node_sparse, density, activate_node_sparse, old_unactivate_nodes);
    % If no new stripped data points are generated, the loop is skipped.
    if isempty(new_add_unactivate_nodes)
        if (verbose)
            disp('new_unactivate_nodes are none, breaking');
        end
        break;
    end
    currIterNdLen = length(old_unactivate_nodes);
    % If the change in core points between two iterations is small (<NdChangeThres), the iteration is stopped.
    if (abs(currIterNdLen - prevIterNdLen) <  NdChangeThres  & iters~=1 )
            if (verbose)
                string_disp = ['stopping peeling since difference between remaining nodes and current is: ',num2str((abs(size(current_data,1) - previous_iteration_data_length)))];
                disp(string_disp);
            end
            break;
    end
    
    prevIterNdLen = currIterNdLen;
    
    [link_core_nodes_indices unactivate_points]= link_core(dist_node, node_ind_point_cell, activate_node_sparse, new_add_unactivate_nodes,density);
    link_border_nodes_indices = new_add_unactivate_nodes;
    border_nodes_indices_per_iteration{end+1} = link_border_nodes_indices;
    core_nodes_indices_per_iteration{end+1} = link_core_nodes_indices;
    activate_data_remat(:,unactivate_points) = 0;
    oneiter_single_gain_mat = single_gain_mat.*activate_data_remat;
    oneiter_single_gain_onedim = reshape(oneiter_single_gain_mat,no_objs*2^no_dims, []);
    density = sum(bsxfun(@times, membership_remat,oneiter_single_gain_onedim))';
    densIters(:,end+1) = density;
    
    % When the density of the boundary points changes dramatically, 
    % it proves that the core region has been found and ends the stripping operation
    if(eps > 0)
            % Find the average density of the boundary points stripped away in this iteration
            newAddBorDens = densIters(new_add_unactivate_nodes,end-1);
            meanBorderDens(end+1) = mean(newAddBorDens);
            
            % Iterate at least three times
            if (length(meanBorderDens) > 2)
                ratioDiff = (meanBorderDens(end) / meanBorderDens(end-1)) - (meanBorderDens(end-1) / meanBorderDens(end-2));
                if verbose
                     string_disp = ['mean border ratio difference:',num2str(ratio_diff)];
                     disp(string_disp)
                end

                if (ratioDiff > eps)
                    if verbose
                        disp('mean border ratio is larger than set value, stopping peeling');
                    end
                    break;
                end  
            end
    end
    
    iters = iters+1;
end

coreNd = setdiff(1:no_node_sparse,old_unactivate_nodes);

distLast = (dist_node.*activate_node_sparse)';
distLast(distLast==0) = Inf;

cluster_lists = mergeCoreNd(distLast,coreNd);
cluster_index = 1;
ndClusters = zeros(1,no_node_sparse);
for l = cluster_lists
    l = l{1};
    if length(l) < min_cluster_size
        continue;
    end
    ndClusters(l) = cluster_index;
    cluster_index = cluster_index + 1;
end

for i = length(border_nodes_indices_per_iteration):-1:1
    link_core_indices = core_nodes_indices_per_iteration{i};
    link_border_indices = border_nodes_indices_per_iteration{i};
    ndClusters(link_border_indices) = ndClusters(link_core_indices);
end


% Node to data point mapping
nearest_node =  round(data);
[~,Locb] = ismember(nearest_node,node_sparse,"rows");
memPtNd = bsxfun(@eq,Locb,1:no_node_sparse);

ptClusters = memPtNd*ndClusters';

end

