% fname = '20211119_cla_cpr.mwk2';
% d = MW_readFile(fname, 'include', {'#stimDisplay'}, 'dotPositions');
% Find confidence interval for noise dots
clear all; close all;
%% At DPZ
%% On Mac, Finder, Command K, cifs://dpz.lokal/dpz
% addpath //Volumes/dpz/userinterchange/KNguyen/ 
%% On windows
% addpath Z:\userinterchange\KNguyen
% In case the network fails:
%load('//dpz.lokal/dpz/userinterchange/Felix Schneider/forKate/dots_only.mat')
%load('\\172.17.14.171\userinterchange\Felix Schneider\forKate\dots_only.mat')

addpath /Users/katenguyen/Documents/DPZ/MATLAB/Data/Solo/ % DPZ windows
addpath /Users/katenguyen/Dropbox/DPZ/MATLAB/Data
addpath /Users/katenguyen/Dropbox/DPZ/MATLAB/Programs
addpath /Users/KateNguyen/Documents/Dropbox/DPZ/MATLAB/Programs
addpath /Users/katenguyen/Documents/DPZ/solo-CPR-model/programs

load('20220128_nak_CPRsolo_block1_fxs.mat')
%load('20220128_nak_CPRsolo_block2_fxs.mat')
%load('forKate.mat')

idx.dot         = d.event == 'STIM_RDP_dot_positions';
dp_all          = d.value(idx.dot);
dp_all_ts       = d.time(idx.dot);

idx.coh         = d.event == 'STIM_RDP_coherence';
dc              = d.value(idx.coh);
dc_ts           = d.time(idx.coh);

idx.rdp_dir     = d.event == 'STIM_RDP_direction';
rdp_dir         = d.value(idx.rdp_dir);
d_rdp_dir_ts    = d.time(idx.rdp_dir); % frame time stamp

dc_mat = cell2mat(dc); % convert cell of coherence to matrix of coherence
dc_unique = unique(dc_mat);

for i = 1:length(dc_unique)
    bl_unique = find(dc_mat == dc_unique(i));
    [a_unique{i},a_stim{i},coh_val] = scprm_programs_extract_func(bl_unique,dc,dc_mat,dc_ts,rdp_dir,d_rdp_dir_ts,dp_all,dp_all_ts);
    coh_cell{i} = coh_val*ones(1,length(a_stim{i}));
end

for i = 1:length(a_stim)
    a_diff{i} = (circ_dist(a_stim{i},a_unique{i}));
end

a_diff_mat = cell2mat(a_diff);

a_diff_nan = a_diff;
a_unique_nan = a_unique;
for i = 1:length(a_diff)
    i_Nan = find(isnan(a_diff{i}) == 1);
    a_diff_nan{i}(i_Nan) = [];
    a_unique_nan{i}(i_Nan) = [];
end

for i = 1:length(a_unique)
    a_mean(i) = mean(a_diff{i},'omitnan');
    a_std(i) = std(a_diff{i},'omitnan');
    coh_mat(i) = coh_cell{i}(1);
end


figure
errorbar(coh_mat,rad2deg(a_mean),rad2deg(a_std),'LineWidth',3)
xlabel('Coherence')
ylabel('Stimulus - Veridical directions (in degrees)')
xlim([0 0.85])
set(gca,'FontSize',30)