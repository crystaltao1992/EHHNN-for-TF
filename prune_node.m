function [B, BB, stem_B, adjacency_matrix, id_layer, id_var_bb, weights, lof, err, stds]=prune_node(B, stem_B, id_layer, id_var_bb, x, y, lambda, parameters)
% �����֦���� 
%

%% parameters
%  complexity penalty
penalty = parameters.penalty;
% lasso parameters
gamma = parameters.gamma;
rho = parameters.rho;
precision = parameters.precision;
quiet = parameters.quiet;

%% ������������Ӿ��� ��ij��Ԫ�ر�ʾ������i->j
num_nodes = size(stem_B,1);

num_connection = nnz(stem_B);
row_indices = zeros(num_connection, 1);
col_indices = zeros(num_connection, 1);
counter = 1;
for kk=1:num_nodes
    for jj = 1:2
        vkk=stem_B(kk,jj);
        if vkk > 0
            row_indices(counter) =  vkk;
            col_indices(counter) =  kk;
            counter = counter + 1;
        end
    end
end
adjacency_matrix = sparse(row_indices, col_indices, 1, num_nodes, num_nodes+1);
%% ��ѵ�����ݣ�����ÿ���ڵ��ֵ��
node_values = cal_node_value(B, stem_B, x);

%% ����ÿ���ڵ㵽������Ȩ��
lambda = lambda * sqrt(2*log10(num_nodes + 1));   % �����������нڵ��������
weights = lasso(node_values, y, lambda, rho, gamma, quiet); % ������ѵ��ÿ���ڵ㵽����ڵ��Ȩ��
weights_of_constant = weights(1);
weights_of_nodes = weights(2:end);

%% ɾ���ڵ�
if sum(abs(weights_of_nodes)) > precision
    index_active_node = abs(weights_of_nodes) > precision;
    adjacency_matrix(index_active_node, num_nodes + 1) = 1;   % the last column
    % ��ɾ������==0�Ľڵ�
    out_fan = sum(adjacency_matrix, 2);
    rem_index = find(out_fan > 0);  % Ҫ�����Ľڵ���
    
    num_node_remained = length(rem_index);
    while num_node_remained < num_nodes
        num_nodes = num_node_remained;
        adjacency_matrix = adjacency_matrix(rem_index, [rem_index',end]);
        % �����������
        B = B(rem_index);
        id_layer = id_layer(rem_index);
        id_var_bb = id_var_bb(rem_index);
        weights_of_nodes = weights_of_nodes(rem_index);
        stem_B = zeros(num_node_remained, 2);
        for nn = 1:num_node_remained
            tmp_id = find(adjacency_matrix(:,nn) > 0)';
            if ~isempty(tmp_id)
                if length(tmp_id) ~= 2
                    fprintf('���Ӿ�������������ڵ�����Ϊ2\n')
                    quit()
                end
                stem_B(nn, :) = tmp_id;
            end
        end
        
        out_fan = sum(adjacency_matrix, 2);
        rem_index = find(out_fan > 0);  % Ҫ�����Ľڵ���
        num_node_remained=length(rem_index);
    end
    % ����ȷ���ڵ����ֵ
    node_values = cal_node_value(B, stem_B, x);
    weights = [weights_of_constant; weights_of_nodes];
    BB = node_values(:,2:end);
    
    hat_y = node_values * weights;
    err = norm(hat_y - y)^2 / norm(y - mean(y))^2;
    stds =  std(hat_y-y);
    lof = err / ( 1 - ( num_nodes + 2 + penalty * (num_nodes+1) ) / size(x, 1) )^2;
else % ���нڵ㵽������Ȩ�ض����㣬ɾ�����нڵ�
    lof = 10;
    err = norm(y)^2 / norm(y - mean(y))^2;
    stds =  std(y);
    if weights_of_constant == 0
        BB=[];
    else
        BB = node_values(:, 2:end);
    end
end

