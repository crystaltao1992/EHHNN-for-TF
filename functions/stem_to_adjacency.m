function [ adjacency ] = stem_to_adjacency( stem )
%STEM_TO_ADJACENCY 
%   stem nx2 matrix
    if isempty(stem) || ~ismatrix(stem) || size(stem, 2) > 2
        fprintf('Invalid input!\n')
        return
    end
    
    num_nodes = size(stem, 1);
    
    [col, ~, value] = find(stem);
    
    adjacency = sparse(value, col, ones(size(col)), num_nodes, num_nodes );


