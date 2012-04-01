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


clear all;
close all;
file_pre_name = 'scissors';

mask_color =  [255, 0, 0; 0,   255, 0];                % colors used to label the foreground and back ground
num_cluster = 32;                                               % number of the clusters 
b_save_result = 1;                                              % if save the result

%file_name

file_save_spline = [file_pre_name, '_spline.bmp']; 

file_kmeans = ['data_image_kmeans_', file_pre_name, '.mat'];

image_file_label = [file_pre_name,   '_label.bmp'];         % you must name the file of the user specified image as  "XX_label.bmp"
image_file          = [file_pre_name,    '.bmp'];

% read the labeled image
A_label = imread(image_file_label); 

figure(1)
imshow(A_label, []);


% read the source color image
Img = imread(image_file);

%===============================  parameters   =====================================
mask_color = [255, 0,   0;  0,   255, 0];   %two


%===============================   extract labels  =================================

[height, width] = size(Img(:,:,1));
N = height * width;
Inds = extract_label_from_image(A_label, mask_color);


%===============================   prepare data  ======================================
X = cvrt_img_to_matrix(Img);

x_coord = ones(height, 1) * ( (1: width) ./ double(width) );            % add the spatial information
y_coord = ( (1: height) ./  double(height) )' * ones(1,  width);

x_coord = x_coord';
y_coord = y_coord';

x_c = x_coord(:);                                                       % convert it to a vector
y_c = y_coord(:);                                                       % convert it to a vector 

%===============================   extract them   ======================================

X_all = [X; x_c'; y_c'];                                            % all the pixel are tested, each column is a data point
X_all = X_all';                                                         % each row is a data point!!!!!!!!!!!


%=========================================================================
%clustering the data points;
%=========================================================================

% record the 
K = size(mask_color, 1);                     % the classes labeled
dim = size(X_all, 2);
XKmeans = zeros(K * num_cluster, dim);

%=========================================================================
% %%% clustering the labelled data,  you can use them by uncommenting them  !!!!!
% for i = 1: K
%     X_temp = X_all(Inds{i, 1}, :);                                     % % % obtain the labeled data points in each class 
%     [I, clusters] = kmeans(X_temp, num_cluster);               
%     XKmeans( ((i - 1) * num_cluster + 1) : (i * num_cluster), :) = clusters;
% end
%  %%%%%%% save(file_kmeans, 'XKmeans');
%=========================================================================

% % Here I give you the clusters.   These culsters are generated by the above sentences in the for-loop. 
%  % If you do not want to use the culsters I provide. you can simply  uncomment the sentences 69-73, and of course comment the following two sentences  
load data_image_kmeans_scissors;
num_cluster = floor((size(XKmeans, 1) ) / 2);

%=========================================================================
% construct the regression values
%=========================================================================

Y = ones(K * num_cluster, 1);
object_value = [-1, 1];    % only work for  K = 2;
for i = 1: K
    Y( ((i - 1) * num_cluster + 1) : (i * num_cluster), 1) = object_value(i); 
end

%=========================================================================
% Spline  regression 
%=========================================================================
tic;
X_test_new =  spline_regression(XKmeans, Y,  X_all);                                % store row-by-row for different methods
time_here = toc
MM = post_process_image_for_two_class_regression(X_test_new, width, height, 1);

figure(2)
imshow(MM, []); 

if  b_save_result == 1
     imwrite(MM, file_save_spline,'bmp');
end


return;







