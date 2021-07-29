clear
close all

dbstop if error
%% load the data

% load('fried_200.mat');  % x_train,y_train,x_test,y_test
% load('fried0113.mat');
load('fried_5000.mat')
config_file = 'config.ini';

% miny = min(A(:,end));
% maxy = max(A(:,end));
% A(:,end) = (A(:,end)-miny)/(maxy-miny);  %you can finish this step by modifyling the data file


%% Parameter Initialiation
parameters = init_par(config_file);
penalty = parameters.penalty;  % complexity penalty
num_train = parameters.num_train;  % number of training
percent = parameters.percent; % percentage of training data
lambda = parameters.lambda;
%% Training set and Test set
% [x_train, y_train, x_validate, y_validate, x_test, y_test] = validate_data(A, percent);

dim = size(x_train,2);
dim_y = size(y_train,2);
num_sample=size(x_train,1);
%%


adjacency_matrices = cell(num_train, 1);
stem_BBs = cell(num_train, 1);
Bs = cell(num_train, 1);
Bs_first = cell(num_train, 1);
layer_indices = cell(num_train, 1);
weights_all = cell(num_train, 1);
lofs = zeros(num_train, 1);
errs = zeros(num_train, 1);
stds_all = zeros(num_train, 1);
f_ehh_tt=zeros(num_sample,num_train);
parameters.structure=0;%20;
% restriction{1}=[1 2 3 4 5];
% restriction{2}=[1 2];
indices = crossvalind('Kfold',num_sample,num_train);
for TT = 1:num_train ;%num_train   %num_train is the training times
%     if TT<floor(num_train/2)+1
%         parameters.structure=0;
%     else
%         parameters.structure=20;
%     end
    id_leave = (indices == TT); 
    id_left = ~id_leave;
    x=x_train(id_left,:);
    y=y_train(id_left);
    [B, weights, id_var_bb, stem_B, adjacency_matrix, id_layer, lof, err, stds, lambda_opt] = forward(x, y, parameters);
%     [B, weights, id_var_bb, stem_B, adjacency_matrix, id_layer, lof, err, stds, lambda_opt]=forward_restricted(x_train, y_train, parameters, restriction);
    
    
    num_nodes=size(stem_B,1);
    pos_row_id = find(stem_B(:,1)>0);  %positive row index, the rows for the first hidden layer are zero
    if isempty(pos_row_id)   % all the neurons are in the first hidden layer
        num1layer = num_nodes;
    else
        num1layer = num_nodes - length(pos_row_id);  % number of nodes in the first hidden layer
    end
    B_first = cell2mat(B(1:num1layer));  % basis function matrix in the first hidden layer
    
    
    
    Bs_first{TT} = B_first;
    adjacency_matrices{TT} = adjacency_matrix;
    stem_BBs{TT} = stem_B;
    Bs{TT} = B;
    layer_indices{TT} = id_layer;
    LAYERS{TT}=[B, num2cell(id_layer,2), num2cell(stem_B,2)];
    LAYERS=LAYERS';
    weights_all{TT} = weights;
    
    lofs(TT)  = lof;
    errs(TT) = err;
    stds_all(TT) = stds;
    f_ehh_TT(:,TT)=cal_node_value(B,stem_B,x_train)*weights;
    
    y_test_pre=cal_node_value(B,stem_B,x_test)*weights;
    err_test(TT)=sqrt(1/length(y_test)*norm(y_test_pre-y_test)^2);%/norm(y_test-mean(y_test))^2;
end

P=2*f_ehh_TT'*f_ehh_TT;
P=(P+P')/2;
q=-2*f_ehh_TT'*y_train;
r=y_train'*y_train;
lb=zeros(num_train,1);
ub=1e6*zeros(num_train,1);
[ratio,~]=quadprog(P, q, [], [],[],[],lb, []);

% ratio = ones(num_train, 1) / num_train;
id_non0=find(ratio>1e-6);
[layers, weights] = merge_net2(LAYERS(id_non0,:), weights_all(id_non0,:), ratio(id_non0), parameters);
B_tt=layers(:,1);
stem_B_tt=cell2mat(layers(:,3));
id_layer_tt=cell2mat(layers(:,2));
% id_var_bb_tt=B_tt(:,2);
yt=cal_node_value(B_tt, stem_B_tt, x_test)*weights;
rmse_test=sqrt(1/length(y_test)*norm(yt-y_test)^2)

[sigma, minusgcv] = anova_ehh(layers, weights, x_train, y_train, parameters);
[~,bb]=sort(sigma(:,2),'descend');
[sigma(bb,:),minusgcv(bb,2)]

%--plot, is available when resctricted forward is run----
var_indices=0:0.05:1;
x_ori=[ones(length(var_indices),5)*0.5];
xlabel_set={'x1','x2','x3','x4','x5'};
close all
for pt=1:5
    x_plot=x_ori;
    x_plot(:,pt)=var_indices';
    y_plot=cal_node_value(B_tt, stem_B_tt, x_plot)*weights;
%     x_plot=(x_plot-0.5)*2;
    y_ori = 10*sin(pi*x_plot(:,1).*x_plot(:,2))+20*(x_plot(:,3)-0.5).^2+10*x_plot(:,4)+5*x_plot(:,5);
    figure
    plot(var_indices,y_plot,'g','linewidth',1.5)
    xlabel(xlabel_set{pt},'fontsize',14)
    hold on
    plot(var_indices,y_ori,'r:','linewidth',1.5)
end
[x1,x2]=meshgrid(var_indices);
length_y=length(x1(:));
x_plot=[x1(:),x2(:),repmat(0.5, length_y, 3)];
y_plot=cal_node_value(B_tt, stem_B_tt, x_plot)*weights;
y_ori = 10*sin(pi*x_plot(:,1).*x_plot(:,2))+20*(x_plot(:,3)-0.5).^2+10*x_plot(:,4)+5*x_plot(:,5);
y_plot=reshape(y_plot,size(x1,1),size(x1,1));
y_ori = reshape(y_ori,size(x1,1),size(x1,1));
% colormap([1 0 0;0 0 1]) %red and blue
figure
surf(x1,x2,y_plot,'edgecolor','g')
hold on
surf(x1,x2,y_ori,'edgecolor','r')
xlabel('x1','fontsize',14)
ylabel('x2','fontsize',14)

