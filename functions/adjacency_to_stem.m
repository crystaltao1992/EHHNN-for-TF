function [ stem ] = adjacency_to_stem( adjacency )
%UNTITLED2 此处显示有关此函数的摘要
%   此处显示详细说明

if isempty(adjacency) 
    stem = [];
    return
end

num_nodes = size(adjacency, 1);
stem = zeros(num_nodes, 2);
for i = 1:num_nodes
    previous = find(adjacency(:, i));
    if isempty(previous) 
        stem(i, :) = [0 0];
    elseif length(previous) ~= 2
        fprintf('one node has not exactly two inputs\n')
        stem = [];
        return
    else
        stem(i, :) = previous(:)';
    end
end

