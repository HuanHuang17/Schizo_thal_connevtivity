%%% thalamus gradient
%%% Leo, 2024/07/05
clear,clc;

addpath(genpath('~\functions'));

% maskdir
maskdir = '~\mask';
% inputDir
inputSubDir = '~\FunImgARWC'; 
% outputDir
outputDir = '~\out'; 
if ~isdir(outputDir)
    mkdir(outputDir);
end

% Gray mask 
CortexMaskFile = [maskdir, filesep, 'cortex_2mm.nii'];
% ROI mask
thaMaskFile = [maskdir, filesep, 'thalamus_2mm.nii'];
%% inputatlas
[atlas_data, header] = y_Read([maskdir, filesep, 'BN_Atlas_246_2mm.nii']);

atlas_data(atlas_data>210) = 0;
uni1 = unique(atlas_data); uni1(1) = [];
ind_atlas = find(atlas_data~=0);
atlas_data1 = atlas_data(atlas_data~=0);

% Gray mask
[~,cMask]=read(CortexMaskFile); ind_GrayMask=find(cMask);
%ROI mask
[~,thaMask]=read(thaMaskFile); ind_ThaMask=find(thaMask);
% Voxel size in mm
voxelsize=2; 
% % Gaussian smoothness in mm
FWHM=4; 

subjectDir = dir(inputSubDir); subjectDir(1:2) = [];

subID = {};
for i = 1:length(subjectDir)
    subID{i} = subjectDir(i).name;
end
save([outputDir, filesep, 'subID.mat'], 'subID');

%% step1: data preprocess and calc thacortical FC
[fc_all, fc_all_GM] = calc_Thacor_FC(inputSubDir, ind_ThaMask, uni1, ind_GrayMask, FWHM, voxelsize, atlas_data1, outputDir);
fc_all(isnan(fc_all)) = 0;
fc_all_mean = mean(fc_all, 3);
save([outputDir, filesep, 'fc_all_atlas.mat'], 'fc_all', '-v7.3');
save([outputDir, filesep, 'fc_all_mean_atlas.mat'], 'fc_all_mean', '-v7.3');
save([outputDir, filesep, 'fc_all_GM.mat'], 'fc_all_GM', '-v7.3');

%% step2: calc CMAP
calc_CMAP(fc_all, fc_all_mean, inputSubDir, thaMaskFile, outputDir);

%% step3: calc COMAP
inputpreprocess = [outputDir, filesep, 'Preprocessd'];
subFGDir = [outputDir, filesep, 'CMAP'];
outputDirCOMAP = [outputDir, filesep, 'COMAP'];
calc_COMAP(inputSubDir, inputpreprocess, subFGDir, outputDirCOMAP, ind_ThaMask, ind_GrayMask, cMask, header);







