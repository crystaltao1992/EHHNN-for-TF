function [node_values, affine_nodes, A_beta, b_beta] = cal_node_affine(B, stem_B, x)
% For each x, calculate the output of each node
% Once the outputs of the nodes in the first hidden layer are fixed, the node outputs of subsequent layers can be calculated
% input: 
%           B         :  cell vector with size num_knots*1, each element is a vector with size 1*3 [1(fixed), subscript of x, beta]
%           stem_B:  matrix with size num_node*2?stem_B(i,i1), stem_B(i,i2)
%           are the indices of nodes to the node i
%           x         :  sample data, num_data*dim, every row is a sample
% output: 
%           node_values: matrix with size num_data*(num_nodes+1)?f(i,j)=f_j(x_i)

num_nodes=size(stem_B,1);
num_data = size(x,1);
sample_set=1:num_data;
% epsilon=1e-8;  % avoid reaching the boundary

if iscell(B)
    pos_row_id=find(stem_B(:,1)>0);  %positive row index, the rows for the first hidden layer are zero
    if isempty(pos_row_id)   % all the neurons are in the first hidden layer
        num1layer=num_bf;
    else
        num1layer=pos_row_id(1)-1;  % number of nodes in the first hidden layer
    end
    B=cell2mat(B(1:num1layer));  % basis function matrix in the first hidden layer
else
    num1layer=size(B,1);
end

node_values(:,1) = ones(num_data,1);  %constant basis
len_beta=num1layer;
A_beta=[];
b_beta=[];
affine_nodes=zeros(num_data, len_beta+1, num_nodes); %recording the affine functions in each node along all the samples

%% the first hidden layer
for i = 1:num1layer  
    S_i=zeros(1, num1layer);
    S_i(i)=1;
    index_x = B(i, 2);
    beta = B(i, 3);
    posi_sample=find( (x(:,index_x) -beta)>=0);
    nega_sample=setdiff( sample_set, posi_sample );
    node_values(posi_sample, i+1)=x(posi_sample, index_x) - beta;
    node_values(nega_sample, i+1)=0;
    affine_nodes(posi_sample,:,i)=[x(posi_sample,index_x), repmat(-S_i, length(posi_sample),1)]; %x(posi_sample,index_x)-beta(i)
%     affine_nodes(nega_sample,:,i)=repmat([0 0], length(nega_sample),1);
    lower_limit=max(x(nega_sample, index_x));
    upper_limit=min(x(posi_sample, index_x));
    if isempty(posi_sample)
        A_beta=[A_beta; -S_i];
        b_beta=[b_beta; -lower_limit];
    elseif isempty(nega_sample)
        A_beta=[A_beta; S_i];
        b_beta=[b_beta; upper_limit];
    else
        A_beta=[A_beta; S_i; -S_i];
        b_beta=[b_beta; upper_limit; -lower_limit];
    end
%     node_values(:, i+1) = max(x(:, index_x) - beta, 0);
end
%% subsequent layers (according to the ascending order of nodes)
for i = num1layer+1:num_nodes  
    input_1 = stem_B(i, 1);
    input_2 = stem_B(i, 2);
    
    min_sample=find( node_values(:,input_1+1)<= node_values(:, input_2+1)); %indices that the first<= the second
    max_sample=setdiff( sample_set, min_sample );
    node_values(min_sample, i+1)=node_values(min_sample, input_1+1);  %taking the minimum
    node_values(max_sample, i+1)=node_values(max_sample,input_2+1);
    affine_nodes(min_sample,:,i)=affine_nodes(min_sample, :, input_1); %does not contain the constant node
    affine_nodes(max_sample,:,i)=affine_nodes(max_sample, :, input_2);
    
    diff_affine=affine_nodes(:,:,input_1)-affine_nodes(:,:,input_2);

    Amin=diff_affine(min_sample, 2:end);  %input_1<=input_2 Amin*beta<=bmin
    bmin=-diff_affine(min_sample,1);
    Amax=-diff_affine(max_sample, 2:end);%input_1>=input_2 Amax*beta>=bmax
    bmax=diff_affine(max_sample,1);
    A_beta=[A_beta;Amin;Amax];
    b_beta=[b_beta;bmin;bmax];
%     node_values(:,i+1) = min(node_values(:,input_1+1), node_values(:,input_2+1));
end
    
