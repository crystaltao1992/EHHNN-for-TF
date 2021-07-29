function  sigma_processed  = process_sigma_tao( sigma, num_setting )
% to conduct spatial-temporary analysis, by qinghua tao
num_lag = num_setting.num_lag;
num_attr = num_setting.num_attr;
num_station = num_setting.num_station;
% initialization
attr = zeros(1, num_attr);
sigma_lag = zeros(1, num_lag);
station = zeros(1, num_station);


l_i = size(sigma, 1);

for i = 1:l_i
    num = sigma(i,1);
    sigma_value = sigma(i,2);
    temp_lag = mod(num,num_lag);
    if temp_lag ~= 0
        sigma_lag(temp_lag) = sigma_lag(temp_lag) + sigma_value;
    else
        sigma_lag(num_lag) = sigma_lag(num_lag) + sigma_value;
    end
    temp_station = ceil(num/30);
    station(temp_station) = station(temp_station) + sigma_value;
%     if temp_floor ~= 0
%         station(temp_station) = station(temp_station) + sigma_value;
%     else
%         station(temp_station) = station(temp_station) + sigma_value;
%     end
    num = mod(num,30);
    if num ~= 0
        temp_attr = ceil(num/10);
        attr(temp_attr) = attr(temp_attr) + sigma_value;
    else 
        attr(num_attr) = attr(num_attr) + sigma_value;
    end
%     if temp_attr ~= 0
%         attr(temp_attr) = attr(temp_attr) + sigma_value;
%     else 
%         attr(3) = attr(3) + sigma_value;
%     end
    
end
sigma_processed.sigma_lag = sigma_lag;
sigma_processed.sigma_attr = attr;
sigma_processed.sigma_station = station;

% def output_list(ls):
%     ls = ['%.4f'%each for each in ls]
%     res = '|'.join(ls)
%     return res
% 
% with open('sigma_info','r') as f:
%     attr=[0,0,0]
%     sigma_lag=[0 for i in range(10)]
%     station=[0 for i in range(7)]
%     for i,line in enumerate(f):
%         num,sigma = line.strip().split()
%         num,sigma = int(num),float(sigma)
%         num-=1
%         sigma_lag[num%10] += sigma
%         station[num//30] +=sigma
%         num = num%30
%         attr[num//10]+=sigma
%     print('Attr:',output_list(attr))
%     print('lag:',output_list(sigma_lag[::-1]))
%     print('station:',output_list(station))


end

