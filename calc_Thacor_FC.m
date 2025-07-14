function [c_all, c_all_GM] = calc_Thacor_FC(inputSubDir, ind_ThaMask, uni1, ind_GrayMask, FWHM, voxelsize, atlas_data1, outputDir)

    subjectDir = dir(inputSubDir); subjectDir(1:2) = [];
    for ii = 1:length(subjectDir)  
        
        c = zeros(length(ind_ThaMask), length(uni1));
        
        [DataLR, ~, ~, ~] = y_ReadAll([inputSubDir, filesep, subjectDir(ii).name]);
        TLR=size(DataLR,4); % Number of time points
    
        Data_Tha = zeros(TLR,length(ind_ThaMask));
        Data_Gray = zeros(TLR,length(ind_GrayMask));
    
        for i=1:TLR      
            DataLR(:,:,:,i)=imgaussfilt3(DataLR(:,:,:,i),FWHM/voxelsize/2.355);
            tmp=DataLR(:,:,:,i);
            Data_Tha(i,:)=tmp(ind_ThaMask);
            Data_Gray(i,:)=tmp(ind_GrayMask);  
        end
    
        %% delete the nan voxel
        Data_Gray(isnan(Data_Gray)) = 0;
        Data_Tha(isnan(Data_Tha)) = 0;
        T = TLR;
    
        % Perform Wishart filter. Glasser et al. 2016
        fprintf('Wishart filtering\n')
        DEMDT=1; %Use 1 if demeaning and detrending (e.g. a timeseries) or -1 if not doing this (e.g. a PCA series)
        VN=1; %Initial variance normalization dimensionality
        Iterate=2; %Iterate to convergence of dim estimate and variance normalization
        NDist=2; %Number of Wishart Filters to apply (for most single subject CIFTI grayordinates data 2 works well)  
    
        OutCC=icaDim(Data_Tha',DEMDT,VN,Iterate,NDist);
        Data_Tha=OutCC.data';
        OutGray=icaDim(Data_Gray',DEMDT,VN,Iterate,NDist);
        Data_Gray=OutGray.data';
            
        % Demean and std
        Data_Tha=detrend(Data_Tha,'constant');Data_Tha=Data_Tha./repmat(std(Data_Tha),T,1); %remove mean and make std=1
        Data_Gray=detrend(Data_Gray,'constant');Data_Gray=Data_Gray./repmat(std(Data_Gray),T,1); %remove mean and make std=1
        
        outData = [outputDir, filesep, 'Preprocessd', filesep, subjectDir(ii).name];
        if ~isdir(outData)
            mkdir(outData);
        end
        save([outData, filesep, 'Data_Tha.mat'], 'Data_Tha');
        save([outData, filesep, 'Data_Gray.mat'], 'Data_Gray');

        %%  calculate the FC matrix   
        data11_GM = [];
        for j = 1:length(uni1)  
            ind2 = find(atlas_data1==uni1(j));
            Data_Gray_ind = Data_Gray(:, ind2);
            data11 = sum(Data_Gray_ind, 2) ./ length(Data_Gray_ind(1, :)~=0);
            c(:, j) = corr(Data_Tha, data11);   
            data11_GM(:,j) = data11;
        end
        c_all(:, :, ii) = c;
        c_all_GM(:, :, ii) = corr(data11_GM);
    
    end
end