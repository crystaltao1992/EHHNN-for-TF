function [B, BB, id_var_bb, stem_B, id_layer]=ini_basis(x,shares)
% inputs:
%           x: the sample x used for the forward procedure, with size num_data x dim
%           shares: the number of points interpolate in the interval
% outputs:
%           B: (dim + num_knot)*1 cell vector, every element is a vector 1*3, i.e.,[1(fixed value), subscript of x, beta]
%           id_var_bb: (dim + num_knot)vector, each element is the subscript of x
%           BB: matrix with size (dim + num_knot) *num_data, the i-th row is max(0, x in id_var_bb{i} multiplied by beta)
%           stem_B: 
%           id_layer: the index of layer in which the basis functions located
%           Beta: matrix with size num_knot*2, including \beta_{km}, the second column is all \beta_{km}, the first column is the subscript of x
%                            is the candidate bias matrix, dynamically changed


%----------- data preprocessing----------
[num_data, dim] = size(x);
beta_candi = linspace(1/shares, 1-1/shares , shares-1);%Equipartition point in x_i

Beta = [];
BB = x';

tmp_id_var_bb = 1:dim;
% basis
for i=1:dim
    % choose beta_i in x_i
    all_samples_i = x(:,i); % all samples in x_i
    knot=[];
    num_0=0;
    kn0=0;
%     for ii = 1:length(beta_candi)
%         num_ii = length(find(all_samples_i<=beta_candi(ii)));
%         samples_num = num_ii - num_0;
%         if samples_num >= 0.2*num_data  % if 20% of the sample points are in the interval, add a middle point
%             knot = [knot;(kn0+beta_candi(ii))/2;beta_candi(ii)];
%         elseif samples_num >= 0.1*num_data
%             knot = [knot;beta_candi(ii)];
%         end
%         num_0 = num_ii;
%         kn0 = beta_candi(ii);
%     end
    knot=beta_candi';
    Beta = [Beta;  i*ones(length(knot),1) knot];
    tmp_BB = max(repmat(x(:,i), 1, length(knot)) - repmat(knot', num_data, 1), 0)';
    BB = [BB; tmp_BB];
    tmp_id_var_bb = [tmp_id_var_bb  i*ones(1, length(knot))];
end

id_var_bb = num2cell(tmp_id_var_bb');
affines = [ones(dim,1) [1:dim]' zeros(dim, 1); 
               ones( size(Beta,1), 1), Beta];
B = num2cell(affines, 2);

l_bf = length(id_var_bb);
stem_B = zeros(l_bf, 2);
id_layer = ones(l_bf,1);