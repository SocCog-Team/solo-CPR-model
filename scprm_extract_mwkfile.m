% Extract mwk file (dot positions, joystick direction, etc.)

clear all; close all;

addpath /Users/katenguyen/Dropbox/DPZ/MATLAB/Data % directory to the data file

fname = '/Users/katenguyen/Documents/DPZ/MATLAB/Data/Solo/20220329_gaa_CPRsolo_block2_fxs.mwk2';

% extract joystick direction, joystick strength and dot position data
d = MW_readFile(fname, 'include', {'IO_joystickDirection','IO_joystickStrength','#stimDisplay'}, 'dotPositions');

% save to the directory
save('/Users/katenguyen/Documents/DPZ/MATLAB/Data/Solo/20220329_gaa_CPRsolo_block2_fxs.mat','d','-v7.3')