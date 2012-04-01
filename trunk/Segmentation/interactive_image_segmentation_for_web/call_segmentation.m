%   version 2.0 --Jun/2010 
%   version 1.0 --Feb/2008 

%   written by Shiming Xiang ( smxiang@gmail.com ) 



function call_segmentation(source_file_path, file_pre_name, save_path, mask_color, num_cluster, b_save_result)

% source_file_path:                  the folder of the source image to be segmented

% file_pre_name:                      the name of the source image, for example, for image '012.bmp',   file_pre_name = '012'; 
%                                             Note that, corresponding to '012.bmp', you must provide '012_label.bmp'  recording  the user labeled pixels.

% save_path:                            the folder for saving the results
%                                             default:   = source_file_path

% mask_color:                         for foreground/background segmetnation,  you can set this parameter as:      mask_color = [255, 0,   0; 0,   255, 0]; 
%                                            default:    = [255, 0, 0; 0,   255, 0];

% num_cluster:                        number of clusters
%                                            default:    = 30;

% b_save_result:                     =1, the segmentation will be saved;   NOT saved, otherwise.
%                                            default:   = 1;

%  ============================================================
% a demo ---- How to use
%  ============================================================

% clear all;
% close all;
% file_pre_name = '012';
% source_file_path = 'E:\xsmmatlab\segmentation_spline\img_seg_results\temp_2010\';
% save_path =  'E:\xsmmatlab\segmentation_spline\img_seg_results\temp_2010\';
% mask_color =  [255, 0, 0; 0,   255, 0];
% num_cluster = 30;
% b_write_result = 1;
% call_segmentation(source_file_path, file_pre_name, save_path, mask_color, num_cluster, b_write_result);

 
if nargin < 3 
    save_path = source_file_path;
    mask_color =  [255, 0, 0; 0,   255, 0];
    num_cluster = 30;
    b_write_result = 1;
elseif nargin < 4 
    mask_color =  [255, 0, 0; 0,   255, 0];
    num_cluster = 30;
    b_write_result = 1;
elseif nargin < 5
    num_cluster = 30;
    b_write_result = 1;
elseif nargin < 6
    b_write_result = 1;
end

% source file
image_file          = [source_file_path,   file_pre_name,    '.bmp'];
image_file_label = [source_file_path,   file_pre_name,   '_label.bmp'];         % you must name the file of the user specified image as  "XX_label.bmp"

%  to save the results
file_save_spline = [save_path, file_pre_name, '_spline.bmp']; 


% read the labeled image
A_label = imread(image_file_label); 

figure(1)
imshow(A_label, []);

% read the source color image
Img = imread(image_file);

%=============================   extract labels  =================================

[height, width] = size(Img(:,:,1));
N = height * width;
Inds = extract_label_from_image(A_label, mask_color);


%========================   prepare data  ======================================
X = cvrt_img_to_matrix(Img);

x_coord = ones(height, 1) * ( (1: width) ./ double(width) );            % add the spatial information
y_coord = ( (1: height) ./  double(height) )' * ones(1,  width);

x_coord = x_coord';
y_coord = y_coord';

x_c = x_coord(:);                                                       % convert it to a vector
y_c = y_coord(:);                                                       % convert it to a vector 

%========================   extract them   ======================================

X_all = [X; x_c'; y_c'];                                            % all the pixel are tested, each column is a data point
X_all = X_all';                                                         % each row is a data point 


%=========================================================================
%clustering the data points;
%=========================================================================

% record the 
K = size(mask_color, 1);                     % the classes labeled
dim = size(X_all, 2);
XKmeans = zeros(K * num_cluster, dim);

% clustering the labelled data
for i = 1: K
    X_temp = X_all(Inds{i, 1}, :);                           % obtain the labeled data points in each class 
    
    [I, clusters] = kmeans(X_temp, num_cluster);    %  note that K-means may yield empty clusters. This is a problem of K-means algorithm. This may be replaced by other clustering algorithms. 
                                                                             % also note that K-means may generate different clusters  if running " kmeans"  in different time. 
                                                                             % Due to this factor, the segmetnation results may be (sightly) different from each other   !!! !!! !!! !!! !!! 
                                                                             % Of course, it is necessary to provide a  robust K-means algorithm. 
                                                                             
    XKmeans( ((i - 1) * num_cluster + 1) : (i * num_cluster), :) = clusters;
end

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







