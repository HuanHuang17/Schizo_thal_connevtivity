function calc_CMAP(fc_all, fc_mean, inputDirSub, thaMaskFile, outputDir)

    subName1 = dir(inputDirSub); subName1(1:2) = [];
    
    % ROI mask
    [thaMaskFile, header] = y_Read(thaMaskFile);
    header.dt = [16, 0];

    fc_all(isnan(fc_all)) = 0;
    fc_mean(isnan(fc_mean)) = 0;
    fc_all(isinf(fc_all)) = 0;
    fc_mean(isinf(fc_mean)) = 0;
    ind_emp = [];
    
    Gref = GradientMaps();
    Gref = Gref.fit(fc_mean,'sparsity', 90);
    
    for s = 1:size(fc_all, 3)
        gm1 = GradientMaps('alignment', 'pa');
        sub1 = fc_all(:, :, s);
        try
            gm1 = gm1.fit( sub1, 'reference', Gref.gradients{1});
        catch
            gm1 = [];
            gm1.gradients{1} = [];
            ind_emp(s) = 1;
        end
    
        for i = 1:size(gm1.gradients{1}, 2)
    
            data1 = gm1.aligned{1}(:,i);
            tmp1 = zeros(size(thaMaskFile));
    
            tmp1(thaMaskFile~=0) = data1;
            out1 = [outputDir, filesep, 'CMAP', filesep, subName1(s).name];
            if ~isdir(out1)
                mkdir(out1);
            end       
            y_Write(tmp1, header, [out1, filesep, 'CMAP_', num2str(i), '.nii']);
        end
    end

end