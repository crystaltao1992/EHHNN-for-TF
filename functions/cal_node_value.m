function node_values = cal_node_value(B, stem_B, x)
% For each x, calculate the output of each node
% Once the outputs of the nodes in the first hidden layer are fixed, the node outputs of subsequent layers can be calculated
% input: 
%           B         :  cell vector with size num_knots*1, each element is a vector with size 1*3 [1(fixed), subscript of x, beta]
%           B also can be matrix
%           stem_B:  matrix with size num_node*2?stem_B(i,i1), stem_B(i,i2)
%           are the indices of nodes to the node i
%           x         :  sample data, num_data*dim, every row is a sample
% output: 
%           node_values: matrix with size num_data*(num_nodes+1)?f(i,j)=f_j(x_i)

num_nodes=size(stem_B,1);
num_data = size(x,1);

if iscell(B)
    pos_row_id = find(stem_B(:,1)>0);  %positive row index, the rows for the first hidden layer are zero
    if isempty(pos_row_id)   % all the neurons are in the first hidden layer
        num1layer = num_nodes;
    else
        num1layer = num_nodes - length(pos_row_id);  % number of nodes in the first hidden layer
    end
    B = cell2mat(B(1:num1layer));  % basis function matrix in the first hidden layer
else
    num1layer=size(B,1);
end

node_values(:,1) = ones(num_data,1);  %constant basis

%% the first hidden layer
for i = 1:num1layer  
    index_x = B(i, 2);
    beta = B(i, 3);
    node_values(:, i+1) = max(x(:, index_x) - beta, 0);
end
%% subsequent layers
for i = num1layer+1:num_nodes  
    input_1 = stem_B(i, 1);
    input_2 = stem_B(i, 2);
    node_values(:,i+1) = min(node_values(:,input_1+1), node_values(:,input_2+1));
end
    
