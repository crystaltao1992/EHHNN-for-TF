function [x_train, y_train,  x_test, y_test] = get_data(A, percent)
% choose train data, validation data, test data
% input:
%       A                :   sample matrix, every row is a sample
%       percent: vector with [percentage of train data,  percentage of test data]

    n_tol = size(A,1);  %total number of sample
    n_train = floor(n_tol*percent(1));  %number of train data
    n_test=n_tol-n_train;  % number of test data

%     id = randperm(n_tol);% randomized sequence of sample
    id = 1:n_tol;
    id_train = id(1:n_train);
    id_test = id(n_train+1:end);

    data_train = A(id_train,:);
    data_test = A(id_test,:);

    x_train = data_train(:,1:end-1);
    y_train = data_train(:,end);


    x_test = data_test(:,1:end-1);
    y_test = data_test(:,end);
    