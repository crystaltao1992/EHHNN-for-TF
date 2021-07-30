% This the demo for showing ENNHH for TF prediction.
% qtao@esat.kuleuven.be, zueslee.hitsz@foxmail.com
clear;
warning off;
dbstop if error
addpath(genpath(pwd));

%% load the data
data_file= 'traffic_10_1_3_origin.data';% predict T=1
config_file = 'config.ini'; % this can be tuned
A = textread(strcat('./data/' , data_file));  % every row stands for a sample, the last column is y

%% Parameter Initialiation
parameters = init_par(config_file); 
% parameters.structure = [0]; % [0] is the single-layered ehh and can be used for individual variable selections with ANOVA decomposition
penalty = parameters.penalty;  % complexity penalty
num_train = parameters.num_train;  % number of training
percent = parameters.percent; % percentage of training data
lambda = parameters.lambda;
lambda = 0.1; % not too much efforts on tuning can commonly obtain good prediction results for almost cases

%% Dataset splitting 
[x_train, y_train, x_test, y_test] = get_data(A, percent);%validate_data(A, percent);
load('index_variable_selection.mat'); % the demo for the selected variables
x_train = x_train(:, index_variable_selection); 
x_test = x_test(:, index_variable_selection);
x_validate = x_train; y_validate = y_train; % save the original training data, which will be split for generating subnets
x_train_all = x_train; y_train_all = y_train;
n_train = size(x_train_all, 1);

%% Train EHHNN predictor for TF
parameters.lambda = lambda;
adjacency_matrices = cell(num_train, 1);
stem_BBs = cell(num_train, 1);
Bs = cell(num_train, 1);
Bs_first = cell(num_train, 1);
layer_indices = cell(num_train, 1);
weights_all = cell(num_train, 1);
LAYERS = cell(num_train, 1);
lofs = zeros(num_train, 1);
errs = zeros(num_train, 1);
stds_all = zeros(num_train, 1);
yahh = zeros(length(y_test), num_train);
err_test = zeros(num_train, 1);
std_test =  zeros(num_train, 1);

% ------------------------ generating subnets-------------------------
for TT = 1:num_train % generate subnets of EHHNN
    x_train = x_train_all(1:n_train - num_train + TT, :);
    y_train = y_train_all(1:n_train - num_train + TT, :);
    [B, weights, id_var_bb, stem_B, adjacency_matrix, id_layer, lof, err, stds, lambda_opt] ...
        = forward_tao(x_train, y_train, parameters);
    num_nodes=size(stem_B,1);
    pos_row_id = find(stem_B(:,1)>0); 
    if isempty(pos_row_id)  
        num1layer = num_nodes;
    else
        num1layer = num_nodes - length(pos_row_id);  % number of nodes in the first hidden layer
    end
    B_first = cell2mat(B(1:num1layer));  % source neurons
    Bs_first{TT} = B_first;
    adjacency_matrices{TT} = adjacency_matrix;
    stem_BBs{TT} = stem_B;
    Bs{TT} = B;
    layer_indices{TT} = id_layer;
    LAYERS{TT}=[B, num2cell(id_layer,2), num2cell(stem_B,2)];
    weights_all{TT} = weights;

    lofs(TT)  = lof;
    errs(TT) = err;
    stds_all(TT) = stds;

    node_values = cal_node_value(B,stem_B, x_test);
    yahh(:,TT) = node_values*weights;
    err_test(TT) = norm( yahh(:,TT) - y_test )^2 / norm( y_test - mean( y_test ) )^2;
    std_test(TT) = std( yahh(:,TT) - y_test);

    f_ehh_TT(:,TT)=cal_node_value(B,stem_B,x_validate)*weights;
end

% ------------------------ stacking subnets-------------------------
P=2*f_ehh_TT'*f_ehh_TT;
q=-2*f_ehh_TT'*y_validate;
r=y_validate'*y_validate;
lb=zeros(num_train,1);
[ratio,~]=quadprog(P, q, [], [],[],[],lb, []);
% merge the subnets
id_non0=find(ratio>1e-6);
[layers, weights] = merge_net2(LAYERS(id_non0,:), weights_all(id_non0,:), ratio(id_non0), parameters );
B_tt = layers(:,1);
stem_B_tt = cell2mat(layers(:,3));
id_layer_tt = cell2mat(layers(:,2));
layers_all = layers;
weights_all = weights;
B_tt_all = B_tt;
stem_B_tt_all = stem_B_tt;

% testing result
ytt = cal_node_value(B_tt, stem_B_tt, x_test)*weights;
test_n =size(y_test,1);
MAE = sum(abs(ytt-y_test))/test_n;
RMSE = sqrt( sum((ytt-y_test).^2)/test_n );
R2 = 1-norm(ytt-y_test)^2/norm(y_test-mean(y_test))^2;

fprintf('Tesingting RMSE: %.3f, MAE: %.3f, R2: %.3f \n', RMSE, MAE, R2)


%% ------------------- anova decomposition ----------------------
% varied variables analysis can be flexibly done here, according to the needs
[sigma_dim, minusgcv_dim] = anova_ehh(layers, weights, x_train_all, y_train_all, parameters);




