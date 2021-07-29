function [B_new, BB_new, stem_B_new, id_var_bb_new, id_layer_new ] = generate_neurons_rand(candidate_combinations,  num_neurons, B, BB, id_var_bb, id_layer)
% generate_neurons_rand generate new neurons randomly
%

    % remove duplication in the candidates
    tt = unique(candidate_combinations, 'rows');
    if size(tt, 1) ~= size(candidate_combinations, 1)
        disp('disp')
        candidate_combinations = tt;
    end
    % permutate the candidates
    shape=size(candidate_combinations);
    rand_comb = randperm(shape(1) );
    num_neurons = min(length(rand_comb), num_neurons);
    comb_choose = candidate_combinations(rand_comb(1:num_neurons),:);
    comb_choose = sortrows(comb_choose);
    
    B_new = cell(num_neurons, 1);
    BB_new = zeros(size(BB, 1), num_neurons);
    stem_B_new = zeros(num_neurons, 2);
    id_var_bb_new = cell(num_neurons, 1);
    id_layer_new = zeros(num_neurons, 1);
    for k = 1:num_neurons
        n1 = comb_choose(k,1);
        n2 = comb_choose(k,2);
        
        B_new{k} = [B{n1}; B{n2}];
        BB_new(:,k) = min(BB(:,n1),BB(:,n2));
        stem_B_new(k,:) = sort([n1, n2]); % sort the constitute neurons 
        id_var_bb_new{k} = union(id_var_bb{n1}, id_var_bb{n2});
        id_layer_new(k) = id_layer(n1) + id_layer(n2);
    end