function [Bt,stem_Bt]=aggregate_adja(stem_B1,B11,id_layer1,coef1,stem_B2,B12,id_layer2,coef2)
%----aggregation of adjacency matrix----
%----A1,A2 can be matrix or stem_B----

num1_max=max(id_layer1);
num2_max=max(id_layer2);
num_max=max(num1_max,num2_max); %num_max is the number of layer after aggregation

stem_B1=[(1:size(stem_B1,1))',stem_B1,coef1(2:end)];
stem_B2=[(1:size(stem_B2,1))',stem_B2,coef2(2:end)];
B11=[(1:size(B11,1))',B11];
B12=[(1:size(B12,1))',B12];

num1=zeros(1,num_max); % ??????????????
for ii=1:num1_max
    num1(ii)=length(find(id_layer1==ii));
end

num2=zeros(1,num_max); % ??????????????
for ii=1:num2_max
    num2(ii)=length(find(id_layer2==ii));
end

stem_Bt=[];
lt=0;

for nn=1:num_max    %build the aggregated neural network layer by layer
    if nn==1
        id_nn1=1:num1(1);
        id_nn2=1:num2(1);
        vec1=B11;
        vec2=B12;
    else
        id_nn1=sum(num1(1:nn-1))+1:sum(num1(1:nn)); % in the original stem_B1
        id_nn2=sum(num2(1:nn-1))+1:sum(num2(1:nn)); % in the original stem_B2
        vec1=stem_B1(id_nn1,:);
        vec2=stem_B2(id_nn2,:);
    end
    if isempty(id_nn1)    %neural network 1 does not has layer nn, all the 2 networks should have layer 1
        stem_Bt=[stem_Bt;vec2];
    elseif isempty(id_nn2)
        stem_Bt=[stem_Bt;vec1];
    else
        [c,i1,i2]=intersect(vec1(:,2:end),vec2(:,2:end),'rows'); %i2 is the index in vec2
        rem_i2=setdiff(1:length(id_nn2),i2); % rem_i2 is the index in vec2
       
        idfrom2=[vec2(i2,1);vec2(rem_i2',1)];     % change the indices in stem_B2 for the nn-th layer
        idto2=[vec1(i1,1);(lt+length(id_nn1)+1:lt+length(id_nn1)+length(rem_i2))'];
        
        stem_temp=stem_B2;
        for ii=1:length(idfrom2)
            b=find(stem_B2==idfrom2(ii));
            stem_temp(b)=idto2(ii);
        end
        stem_B2=stem_temp;
        
        if nn<num1_max
            idfrom1=max(id_nn1)+1:max(id_nn1)+num1(nn+1);% change the indices in stem_B1 for the nn+1-th layer
            idto1=idfrom1+length(rem_i2);
        end
        
        stem_temp=stem_B1;
        for ii=1:length(idfrom1)
            b=find(stem_B1==idfrom1(ii));
            stem_temp(b)=idto1(ii);
        end
        stem_B1=stem_temp;
        
        
       if nn==1
           Bt=[B11(:,2:end);B12(rem_i2,2:end)];
       end
        stem_Bt=[stem_Bt;stem_B1(id_nn1,:);stem_B2(id_nn2(rem_i2),:)];
        lt=size(stem_Bt,1);
    end
end
        
    