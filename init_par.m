function parameters = init_par(config_file)
    fprintf('��ʼ��ȡ�����ļ���%s\n', config_file)
    parameters = ini2struct(config_file);
    fields = fieldnames(parameters);
    for i = 1:length(fields)
        parameters.(fields{i}) = str2num(parameters.(fields{i})); %#ok<ST2NM>
    end
    fprintf('�����ļ���ȡ���\n')