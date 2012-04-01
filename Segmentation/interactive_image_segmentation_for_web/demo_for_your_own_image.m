%==========================================================================
%    CORRESPONDENCE INFORMATION
%    This code is written by Shiming Xiang
% 
%    Address:    National Laboratory of Pattern Recognition (NLPR), Institute of Automation,
%                     Chinese Academy of Sciences, Beijing, 100190, China, 

%    Email:        smxiang@gmail.com

%    This software can be used freely for research purposes.
%    Published reports of research  using this code (or a modified version) should cite the 
%    article that describes the algorithm: 
%
%    Shiming Xiang, Feiping Nie, Chunxia Zhang and Changshui Zhang.
%    Interactive Natural Image Segmentation via Spline Regression. IEEE Transactions on Image Processing, 
%    Volume 18, Issue 7, Pages1623-1632, 2009

%    Comments and bug reports are welcome.  Email to smxiang@gmail.com. 
%    I would also appreciate hearing about how you used this code, 
%    and the improvements that you have made to it, or translations into other languages.    

%   WORK SETTING:
%    This code has been compiled and tested using matlab    7.0
%==========================================================================

%   version 2.0 --Jun/2010 
%   version 1.0 --Feb/2008 

% Rewritten time:   June 7-9, 2010


%  ============================================================
% a demo:   the image is used in the Fig.1 in the paper 
%  ============================================================

clear all;
close all;
file_pre_name = '002';                       

source_file_path = 'E:\xsmmatlab\segmentation_spline\Code_for_Web\'                                % you must change it, and let it point to your own folder of  image files       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
save_path          = 'E:\xsmmatlab\segmentation_spline\Code_for_Web\'                                % you must change it, and let it point to your own folder of  image files       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

mask_color =  [255, 0, 0; 0,   255, 0];                % colors used to label the foreground and back ground
num_cluster = 30;                                               % number of the clusters 
b_save_result = 1;                                              % if save the result
call_segmentation(source_file_path, file_pre_name, save_path, mask_color, num_cluster, b_save_result);

 return;