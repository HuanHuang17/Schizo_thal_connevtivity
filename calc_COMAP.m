function calc_COMAP(inputSubDir, inputpreprocess, subFGDir, outputDir, ind_ThaMask, ind_GrayMask, cMask, header)

    subjectDir = dir(inputSubDir); subjectDir(1:2) = [];
    for ii = 1:length(subjectDir)  
        
        Data_Tha = importdata([inputpreprocess, filesep, subjectDir(ii).name, filesep, 'Data_Tha.mat']);
        Data_Gray = importdata([inputpreprocess, filesep, subjectDir(ii).name, filesep, 'Data_Gray.mat']);
                
        individual_FG1 = y_ReadAll([subFGDir, filesep, subjectDir(ii).name]);
        individual_FG1 = reshape(individual_FG1, size(individual_FG1, 1)*size(individual_FG1, 2)*size(individual_FG1, 3), size(individual_FG1, 4));
        individual_FG1 = individual_FG1(ind_ThaMask, :);

        %% COMAP        
        one1 = ones(size(individual_FG1, 1), 1);
        for ind = 1:size(individual_FG1, 2)
            fg1 = individual_FG1(:, ind);
            X = [one1, fg1];
            beta1 = pinv(X) * Data_Tha';
            map1_all = Data_Gray' * beta1';
            map1 = map1_all(:, 2);
    %         map1 = (map1 - min(map1)) ./ max(map1);          
            data_emp = zeros(size(cMask));
            data_emp(ind_GrayMask) = map1;
            
            out1 = [outputDir, filesep, subjectDir(ii).name];
            if ~isdir(out1)
                mkdir(out1);
            end
            header1 = header;
            header1.dt = [16, 0];
            y_Write(data_emp, header1, [out1, filesep, 'COMAP_', num2str(ind), '.nii']);     
        end
    end
    
end