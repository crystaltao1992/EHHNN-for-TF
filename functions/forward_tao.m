function [B, weights, id_var_bb, stem_B, adjacency_matrix, id_layer, lof, err, stds, lambda_opt]=...
    forward_tao(x, y, parameters)
% modified by Qinghua Tao to accelerate the search for candidate combination
% input
%       x          ---------------- the sample x used for the forward
%                   procedure, with size N x dim, part of the whole samples
%       y          ---------------- the sample y used for the forward
%                   procedure, with size N x 1, part of the whole samples

% output
%       BBf      --- the basis function evaluated at the points
%       Bf       --- the parameters of the basis function
%       coe      --- the coefficient matrix, dim x 1
%--forward growing of the network-----
%--random strategy---

% Parameter initilization
shares = parameters.shares;  % quantile number of each coordinate, the number of points interpolate in the interval
structure_parameter = parameters.structure;  % the parameters for structure definition
lambda = parameters.lambda;% the parameter for the Lasso regression 

start_time = tic;
%% the first layer
[B0, BB0, id_var_bb0, stem_B0, id_layer0] = ini_basis(x,shares);  % not containing the constant basis

B = B0;
BB = BB0';
stem_B = stem_B0;
id_layer = id_layer0;
id_var_bb = id_var_bb0;

%% forward process
num_layer = size(structure_parameter, 2);  % the first layer is not taking into consideration
for layer_index = 2:num_layer+1  % the neurons are added layerwisely
    num_neurons = structure_parameter(layer_index-1);
    % All possible combinations of neurons
    candidate_combinations = []; %the index of id_layer
    layer_index_1 = 1;
    if num_neurons ~= 0 
        while layer_index_1 < layer_index   %possible combinations yielding layer nl
            layer_index_2 = layer_index-layer_index_1;
            index_in_layer_1 = find(id_layer==layer_index_1);  % in the k1-th layer
            index_in_layer_2 = find(id_layer==layer_index_2);
            [index_in_layer_1, index_in_layer_2] = meshgrid(index_in_layer_1,index_in_layer_2);
            index_in_layer_1 = index_in_layer_1(:);
            index_in_layer_2 = index_in_layer_2(:);
            index_combinations = [];
            for index_combination = 1:length(index_in_layer_1)
                index_x_1 = id_var_bb{index_in_layer_1(index_combination)};
                index_x_2 = id_var_bb{index_in_layer_2(index_combination)};
                if ~isempty( intersect(index_x_1,index_x_2)) % if they have common x_i, they will not be combined
                    continue;
                end
                index_combinations = [index_combinations; index_combination];
            end
            candidate_combinations = [candidate_combinations; sort([index_in_layer_1(index_combinations),index_in_layer_2(index_combinations)], 2)];
            layer_index_1 = layer_index_1+1;
        end
        % Choose suitable combinations from all possible combinations (random procedure)   
        [B_new, BB_new, stem_B_new, id_var_bb_new, id_layer_new ] = generate_neurons_rand(candidate_combinations,  num_neurons, B, BB, id_var_bb, id_layer);
    else
        [B_new, BB_new, stem_B_new, id_var_bb_new, id_layer_new ] = generate_neurons_rand(candidate_combinations,  num_neurons, B, BB, id_var_bb, id_layer);
   
    end
    % -------------------------------- 检查每个节点是否有重复x,
    for iii = 1:size(B_new, 1)
        tmp = B_new{iii};
        if length(unique(tmp(:,2))) ~= size(tmp, 1)
            iii;
        end
    end
    % ---------------------------------check if stem_B_new has duplicated rows
    tmp = unique(stem_B_new, 'rows');
    if size(tmp, 1) ~= size(stem_B_new, 1)
        iii;
    end
    % ---------------------------------
    B = [B; B_new];
    BB = [BB BB_new];
    stem_B = [stem_B; stem_B_new];
    id_var_bb = [id_var_bb; id_var_bb_new];
    id_layer = [id_layer; id_layer_new];
end
 
% ------------------check if stem has duplications
tmp = find(stem_B(:, 1));
tmp = stem_B(tmp, :);
tmp1 = unique(tmp, 'rows');
if size(tmp1, 1) ~= size(tmp, 1)
    TT
end
% ------------------------
    
%% Backward process

lof = inf;
B0=B;
stem_B0=stem_B;
id_layer0=id_layer;
id_var_bb0=id_var_bb;

for k = 1:length(lambda)
    [Bk, BBk, stem_Bk, Adjak, id_layerk, id_var_bbk, coefk, lofk, errk, stdsk] = prune_node(B0, stem_B0, id_layer0, id_var_bb0, x, y,lambda(k), parameters);
    
    % ------------------check if stem has duplications
    tmp = find(stem_Bk(:, 1));
    tmp = stem_Bk(tmp, :);
    tmp1 = unique(tmp, 'rows');
    if size(tmp1, 1) ~= size(tmp, 1)
        k
    end
    % ------------------------
    
    execute_prune = 'X';
    if lofk < lof
        B = Bk;
        BB = BBk;
        stem_B = stem_Bk;
        adjacency_matrix = Adjak;
        id_layer = id_layerk;
        id_var_bb = id_var_bbk;
        weights = coefk;
        lof = lofk;
        execute_prune = 'ok';
        lambda_opt=lambda(k);
    end
    fprintf('lambda: %2.2f, error: %6.4f, lof: %6.4f, std: %6.4f, prune? %s \n', lambda(k), errk, lofk, stdsk, execute_prune);
end

%% The output
node_values = cal_node_value(B,stem_B,x);
hat_y = node_values*weights;
err = norm(hat_y - y)^2/norm(y-mean(y))^2;
stds = std(hat_y - y);
time = toc(start_time);

fprintf('Final results: lambda: %2.2f, error: %6.4f, lof: %6.4f, std: %6.4f, ellapsed time: %f \n', lambda_opt, err, lof, stds, time);
