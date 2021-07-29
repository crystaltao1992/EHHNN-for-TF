function [ interaction ] = adajcency2interaction( adjacency )
%ADAJCENCY2INTERACTION transform adjacency matrix to interaction matrix
% input:  
%           adjacency  : n*n , adjacency matrix
% output:
%           interactioin : n*n , interaction matrix

MAX = 100;

if ~ismatrix(adjacency)
    fprintf('unaccepted input type.\n')
    return
end

num = size(adjacency, 1);
if size(adjacency, 2) ~= num
    fprintf('should be square matrix.\n')
    return
end

mapping = cell(num, 1);
tmp = sum(adjacency, 1);
num_first = sum(tmp == 0);
mapping(1:num_first) = num2cell((1:num_first)', 2); % the first layer nodes which are contained in the i-th node

col = zeros(MAX, 1);
row = col;
count = 1;
for i = num_first+1:num
    prev = adjacency(:, i) > 0; % directly consisting nodes 
    
    consist_first = cell2mat(mapping(prev));
    consist_first = unique(consist_first);
    consist_first = consist_first(:); %  containing nodes in the first layer
    
    mapping{i} = consist_first;
    
    num_containing = length(consist_first);
    col(count: (count + num_containing - 1)) = consist_first;
    row(count: (count + num_containing - 1)) = ones(num_containing, 1) * i;
    
    count = count + num_containing;
end

count = count - 1;
if (count < MAX)
    col = col(1:count);
    row = row(1:count);
end
row = [row; (1:num_first)'];
col = [col; (1:num_first)'];

interaction = sparse(col, row, 1, num, num);
