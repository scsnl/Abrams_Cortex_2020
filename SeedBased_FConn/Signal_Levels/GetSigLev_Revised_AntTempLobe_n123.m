clear all; close all; clc

subs = ReadList('/oak/stanford/groups/menon/projects/daa/2015_ACx_parcellation_HCP/data/subjectlist/Orig_Sample_n123/Behav_Movement_Eligble_HCP_ArtRem_n123.txt');

% Using the same order as the Orig script used to plot SCSNL data
roi_names = ReadList('/oak/stanford/groups/menon/projects/daa/2015_ACx_parcellation_HCP/data/imaging/roi/Revised_ROIs_AntTemporalLobe/ROI_List_18_Bilateral_STC_ResultDirNames.txt');
% This gets rid of right-hemi STC and SPT seeds headers
roi_names([9:16 18]) = [];

% This is for creating headers
roi_names_short = ReadList('/oak/stanford/groups/menon/projects/daa/2015_ACx_parcellation_HCP/data/imaging/roi/Revised_ROIs_AntTemporalLobe/ROI_List_18_Bilat_STC_Seeds_ShortNames.txt');
% This gets rid of right-hemi STC and SPT seeds headers
roi_names_short([9:16 18]) = [];


target_roi_dir = '/oak/stanford/groups/menon/projects/daa/2015_ACx_parcellation_HCP/data/imaging/roi/Revised_ROIs_AntTemporalLobe/NIIs';

get_target_fnames = {'01-4mm_Left_BA44_Peak_Neurosynth_Speech_-46_14_22_-46_14_22.nii',...
    '02-4mm_Left_BA45_Peak_Neurosynth_Speech_-48_24_24_-48_24_24.nii',...
    '03-4mm_Left_BA47_Peak_Neurosynth_Speech_-46_26_-6_-46_26_-6.nii',...
    '04-4mm_Left_aSMG_Peak_Neurosynth_Speech_-50_-38_42_-50_-38_42.nii',...
    '05-4mm_Left_PGa_Peak_Neurosynth_Speech_-30_-64_48_-30_-64_48.nii',...
    '06-4mm_Left_PGp_Peak_Neurosynth_Speech_-48_-68_18_-48_-68_18.nii',...
    '04-4mm_PreCentral_Left_-48_-2_40.nii'};

% This is for creating headers
target_names_short = {'Left_BA44',...
    'Left_BA45',...
    'Left_BA47',...
    'Left_aSMG',...
    'Left_PGa',...
    'Left_PGp',...
    'Left_PreCentral'};

% Check target ROIs
for target_i = 1:length(get_target_fnames)
    
    target_names{target_i,1} = get_target_fnames{target_i};
    
    [target, hdr] = cbiReadNifti(fullfile(target_roi_dir, get_target_fnames{target_i}));
    
    target_mask = find(target > 0);
    
    target_length(target_i) = length(target_mask);
    
end

tic

root_dir = '/oak/stanford/groups/menon/projects/daa/2015_ACx_parcellation_HCP/results/restfmri/participants';

for sub_i = 1:length(subs)
    
    for roi_i = 1:length(roi_names)
                
        hcp_fname = fullfile(root_dir, subs{sub_i}, 'visit1/session1/seedfc/resting_state_1/stats_spm8',...
            roi_names{roi_i}, 'stats', 'con_0001.img');
        
        
        [fc_data, hdr] = cbiReadNifti(hcp_fname);
        
        for target_i = 1:length(get_target_fnames)
          
            target_names{target_i,1} = get_target_fnames{target_i};
            
            [target, hdr] = cbiReadNifti(fullfile(target_roi_dir, get_target_fnames{target_i}));
            
            target_mask = find(target > 0);
            
            Polar_Matrix{roi_i}(sub_i, target_i) = nanmean(fc_data(target_mask));

        end
        
    end
end

toc

no_subs = length(subs);

Monster_Matrix = [zeros(no_subs,1)];
Monster_Matrix_Header = cell(1,1);

target_names_short{8} = 'Empty_Row';

% This loop reorders the data so it can be imported into SPSS for RMANOVAs 
for roi_i = 1:length(roi_names)
    
    temp_matrix = Polar_Matrix{roi_i}(:,:);
    
    next_col = size(Monster_Matrix,2) + 1;
    
    Monster_Matrix(:,next_col:next_col+7) = [temp_matrix zeros(no_subs,1)];
    
    Monster_Matrix_Header(:,next_col:next_col+7) = strcat(roi_names_short{roi_i}, '_', target_names_short);
    
end

Monster_Matrix(:,1) = [];
Monster_Matrix_Header(:,1) = [];

% Identify NaNs in results
[nan_row, nan_col] = find(isnan(Monster_Matrix));
% Identify Zeros in results
[zero_row, zero_col] = find(Monster_Matrix == 0);
% DA Changed to add PreCent 05/18/2018
all_zero_cols = 8:8:72;
A = zero_col';
B = all_zero_cols;
[logical_member] = ismember(A,B);
overlap_ind = find(logical_member == 1);
zero_col(overlap_ind) = [];
zero_row(overlap_ind) = [];

% ===================================================
% Create summary stats for plotting
no_subs = length(subs);
mean_monster = mean(Monster_Matrix);
ste_monster = std(Monster_Matrix) ./ sqrt(no_subs);
 
% DA Changed to add PreCent 05/18/2018
% first_col = 1:7:78;
% last_col = 7:7:84;
first_col = 1:8:89;
last_col = 8:8:96;

% This transform a row vector with all summary stats (ie, the
% "mean_monster" variable into a matrix for plotting
for roi_i = 1:length(roi_names)
    
    % DA Changed to add PreCent 05/18/2018
    mean_mat(roi_i,1:8) = mean_monster(first_col(roi_i):last_col(roi_i));
    ste_mat(roi_i,1:8) = ste_monster(first_col(roi_i):last_col(roi_i));
    
end

row_header_mat = roi_names_short';
col_header_mat = target_names_short;

% Delete this last row of zeros
mean_mat(:,8) = [];
ste_mat(:,8) = [];

results_dir = '/oak/stanford/groups/menon/projects/daa/2015_ACx_parcellation_HCP/scripts/restfmri/seedfc/Extract_Sig_Levels/Revised_ROIs_AntTemporalLobe/Extracted_Sig_Levels';
save_fname = fullfile(results_dir, 'MonsterMatrix_HCP_LH_AntTempLobe_n123_Mean_STE.mat');
save(save_fname, 'mean_mat', 'ste_mat', 'row_header_mat', 'col_header_mat', 'Monster_Matrix', 'Monster_Matrix_Header');




