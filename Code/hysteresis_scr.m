% An example how to use hysteresis.m function for modelling and elimination
% of hysteresis effect from Scintrex CG-3M gravity readings from MATLAB
% Command Window or a script file.

% Please edit the path beneath to match your case.
main_folder='D:\Moji_podaci\RADOVI_clanci\2015_Hysteresis\Hysteresis_Code\Data\Input\';

data_files_1=[{[main_folder '4372\010615K2.DAT']}; {[main_folder '4372\020615K2.DAT']}; {[main_folder '4372\030615K2.DAT']}; {[main_folder '4372\050615K2.DAT']}];
data_files_2=[{[main_folder '4373\010615K1.DAT']}; {[main_folder '4373\020615K1.DAT']}; {[main_folder '4373\030615K1.DAT']}];

% output_hyst=hysteresis(in_files, nor, diagram, file_el)
output_hyst_1=hysteresis(data_files_1, 0, 1, 1);
output_hyst_2=hysteresis(data_files_2, 0, 1, 1);

% Comparison of time shifts for simultaneous gravity measurements of two
% gravimeters. (Please note that for June 5th, data of only one gravimeter
% is supplied.)
[output_hyst_1(3:15,1) output_hyst_2(3:end,1)]