function [layers, weights] = merge_net2(LAYERS, WEIGHTS, ratio, parameters )
%MERGE_NET merge multiple PWL networks
%   syntax:
%           [layers, weights] = merge_net2(LAYERS, WEIGHTS, ratio, parameters )
%   input:
%           LAYERS: n*1 struct, fields are B, stem_B, id_layer, n is the
%                   number of sub nets
%           (note: there is another kind of input, i.e.
%           LAYERS: n*1 cell, each element of which is a cell representing
%                   a sub net, and it is a m*3 cell array, the first column
%                   is B, the second column is id_layer, and the last
%                   column is the stem_B , n is the number of sub nets, m is the number of nodes in each sub nets)
%           WEIGHTS: n*1 cell, each of which is the weights of each node in
%                   sub nets to the output
%           ratio : n*1 vector, weights of each sub net
%           parameters: struct, configurations
%   output: 
%           layers: of the same class as the input LAYERS
%           weights: weighted sum weights of each node to the output in the
%               merged net

%% Check the data
if  ~iscell(WEIGHTS)  || ~isvector(ratio)
    fprintf('input arguments are not of the required class.\n')
    return
end

mode = 'cell';
if isstruct(LAYERS)
    num_nets = length(LAYERS);
    tmp = cell(num_nets, 3);
    for i = 1:num_nets
        tmp{i} = [LAYERS(i).B num2cell(LAYERS(i).id_layer, 2) num2cell(LAYERS(i).stem_B, 2)];
    end
    LAYERS = tmp;
    mode = 'struct';
end

num_nets = length(LAYERS);
if length(WEIGHTS) ~= num_nets || length(ratio) ~= num_nets
    fprintf('input arguments are not of the same dimension.\n')
    return
end

% if abs(sum(ratio) - 1) > parameters.precision
%     fprintf('weights of each net are not summed up to 1.\n')
%     return
% end

%% Data in each sub-network
 
num_nodes_separate = zeros(num_nets, 1); % number of nodes in each sub-network
for i = 1:num_nets
    num_nodes_separate(i) = size(LAYERS{i}, 1);
end
num_nodes = sum(num_nodes_separate);
num_accumulated = cumsum(num_nodes_separate);

%% merge all the networks (raw edition)
adjacency = sparse(num_nets);  % simply combine all adjacency matrices of subnets, (diagnal block)
layers = cell(num_nodes, 3);  % list all nodes of subnets, according to net number
weight_constant = 0;
for i = 1:num_nets
    if i == 1
        range = 1:num_accumulated(i);
        layers(1: num_accumulated(1), :) = LAYERS{i};
    else
        range = (num_accumulated(i-1) + 1):num_accumulated(i);
        layers((num_accumulated(i - 1) + 1): num_accumulated(i), :) = LAYERS{i};
    end
    tmp_stem = cell2mat(LAYERS{i}(:, 3));
    adjacency(range, range) = stem_to_adjacency(tmp_stem);
    adjacency(range, num_nodes + 1) = WEIGHTS{i}(2:end) * ratio(i);
    weight_constant = weight_constant + WEIGHTS{i}(1) * ratio(i);
end

corresponding_ids =  cell2mat(layers(:, 2))  ; % [1-layer index]

%%  rearrange nodes so that the adjacent nodes are of the same layer, except the nodes on the boundary
[~, index_rearranged] = sortrows(corresponding_ids);
adjacency = adjacency(index_rearranged, [index_rearranged; end]);
layers = layers(index_rearranged, :);

num_layers = max(corresponding_ids); % number of layers
num_nodes_layerwise = zeros(num_layers, 1); % number of nodes in each layer
for i = 1:num_layers
    range = corresponding_ids == i; 
    num_nodes_layerwise(i) = sum(range);
end
%% find duplicate nodes in the first layer, and merge
node_info_first_layer = cell2mat(layers(1: num_nodes_layerwise(1), 1));
% dbstop if (size(node_info_first_layer, 1) ~= num_nodes_layerwise(1) || size(node_info_first_layer, 2) ~= 3)

[~, i_unique, i_recover] = unique(node_info_first_layer(:, 2:3), 'rows');
num_unique = length(i_unique);
flag = true;
if num_unique < num_nodes_layerwise(1) % duplication occurs
    index_deleted = setdiff(1:num_nodes_layerwise(1), i_unique);
    layers(index_deleted, :) = []; % delete corresponding nodes in node list
    for i = 1:num_unique % merge successors to one node
        index_same = find(i_recover == i); % it is the index in the original adjacency matrix
        if length(index_same) == 1  % unique, no operation
            continue
        end
        value_sum = sum(adjacency(index_same, :), 1);
        % dbstop if (any(value_sum ~= 1))
        adjacency(min(index_same), :) = value_sum;
    end
    adjacency(index_deleted, :) = [];
    adjacency(:, index_deleted) = [];
else % every node in the first layer is unique, then no duplication occurs in the whole merged nets
%     weights = [weight_constant; adjacency(:, end)];
%     stem = adjacency_to_stem(adjacency(:, 1:end-1));
%     layers(:, 3) = num2cell(stem, 2);
    flag = false; 
end

%% find duplicate nodes, layerwise from the second layer
if flag 
    for i = 2:num_layers
        node_index = find(cell2mat(layers(:, 2)) == i); 
    %     dbstop if (length(node_index) ~= num_nodes_layerwise(i))
        reference = zeros(num_nodes_layerwise(i), 2*i);
        for j = 1:num_nodes_layerwise(i)
            reference(j, :) = [sort(layers{node_index(j), 1}(:, 2));   sort(layers{node_index(j), 1}(:, 3))]';
        end
    %     [~, i_unique, i_recover] = unique(reference, 'rows');
        [~, i_unique, i_recover] = unique(adjacency(:, node_index)', 'rows');
        num_unique = length(i_unique);
        if num_unique < num_nodes_layerwise(i) % duplication occurs
            tmp_index = 1:num_nodes_layerwise(i);
            index_deleted = setdiff(tmp_index, i_unique);
            index_deleted = node_index(index_deleted); % index of to-be-deleted nodes in the adjacency matrix
            layers(index_deleted, :) = []; % delete corresponding nodes in node list
            for j = 1:num_unique % merge successors to one node
                index_same = i_recover == j; % it is the index in the sub adjacency matrix
                index_same = node_index(index_same); % it is the index in the original adjacency matrix
                if length(index_same) == 1  % unique, no operation
                    continue
                end
                adjacency(min(index_same), :) = sum(adjacency(index_same, :), 1);
            end
            adjacency(index_deleted, :) = [];
            adjacency(:, index_deleted) = [];
        else %every node in the i-th layer is unique
            continue
        end
        % dbstop if (size(adjacency, 1) ~= size(layers, 1))
        % dbstop if (size(adjacency, 1) ~= (size(adjacency, 2) - 1))
    end
end
%% prepare output
% id_var_bb
for i = 1:size(layers, 1)
    layers{i, 4} = layers{i, 1}(:, 2);
end
% others
weights = [weight_constant; adjacency(:, end)];
if length(weights) ~= (size(layers, 1) + 1)
    disp('error')
end
stem = adjacency_to_stem(adjacency(:, 1:end-1));
if strcmp(mode, 'cell')
    layers(:, 3) = num2cell(stem, 2);
else
    B = layers(:, 1);
    id_layer = cell2mat(layers(:, 2));
    id_var_bb = layers(:, 4);
    layers = struct('B', B, 'id_layer', id_layer, 'stem_B', stem, id_var_bb, 'id_var_bb');
end