function [sigma, minusgcv] = anova_ehh(layers, weights, x, y, parameters)

B = layers(:,1);
id_layer = cell2mat(layers(:,2));
stem_B = cell2mat(layers(:,3));
id_var_bb = layers(:,4);
weights_bb = weights(2:end);

% adjacency = layers(:,5);
% adja = adjacency(:,1:end-1);
% weights = adjacency(2:end:end);
penalty = parameters.penalty;

node_values = cal_node_value(B, stem_B, x);

f = node_values*weights;


num_neuron = length(weights)-1;



% %--univariate summation---
% id_layer1 = find(id_layer==1);
% id_var_bb1 = cell2mat(id_var_bb(id_layer1));
% id_var1 = unique(id_var_bb1);
% sigma = zeros(length(id_var1),2);
% minusgcv = sigma;
% 
% for ii = 1:length(id_var1)
%     var_ii = id_var1(ii);
%     id = find(id_var_bb1 == id_var1(ii));
%     Bii = B(id);
%     stem_Bii = stem_B(id,:);
%     weights_ii = weights_bb(id);
%     fii = node_values(:,id+1)*weights_ii;
%     num_nodes = num_neuron-length(id);
%     sigma(ii,:) = [var_ii, sqrt(var(fii))];
%     fminus = f-fii;
% %     BBminus = node_values(:, 2:end);
% %     BBminus(:,id)=[];
%     cm = num_nodes+1;%trace(BBminus*inv(BBminus'*BBminus)*BBminus')+1;
%     minusgcv(ii,:) = [var_ii,norm(fminus-y)^2 / ( 1 - ( cm + penalty * num_nodes ) / size(x, 1) )^2/norm(y-mean(y))^2];%size(x,1)
% end
% 
% sigma1 = sigma;
% minusgcv1 = minusgcv;

%---anova decomposition----
num_layer = max(id_layer);
id_neuron = zeros( num_neuron, num_layer );
id_neuron = logical( id_neuron );
for l = 1:num_layer
    id_neuron(:, l) = (id_layer==l);
end
id_var_bb1 = cell2mat( id_var_bb( id_neuron(:,1)));
id_var = unique(id_var_bb1);
[ adjacency ] = stem_to_adjacency( stem_B );
[ interaction ] = adajcency2interaction( adjacency );

sigma = zeros(length(id_var),2);  
minusgcv = sigma;
for ii = 1:length(id_var)
    var_ii = id_var(ii);  % the particular variable var_ii
    id1 = (id_var_bb1 == var_ii);  % the rows var_ii lying in (1st hidden layer)
    vd = interaction(id1, :);
    id_ii = (vd == 1);  %id_ii is the neurons containing var_ii
    id_ii = logical(sum(id_ii, 1));
    Bii = B(id_ii);
    stem_Bii = stem_B(id_ii,:);
    weights_ii = weights_bb(id_ii);
    BB = node_values(:, 2:end);
    fii = BB(:,id_ii)*weights_ii;
    
    num_nodes = num_neuron-sum(id_ii);
    sigma(ii,:) = [var_ii, sqrt(var(fii))];
    fminus = f-fii;
%     BBminus = node_values(:, 2:end);
%     BBminus(:,id)=[];
    cm = num_nodes+1;%trace(BBminus*inv(BBminus'*BBminus)*BBminus')+1;
    minusgcv(ii,:) = [var_ii,norm(fminus-y)^2 / ( 1 - ( cm + penalty * num_nodes ) / size(x, 1) )^2/norm(y-mean(y))^2];%size(x,1)
end
   
% sigma1 - sigma;

    
