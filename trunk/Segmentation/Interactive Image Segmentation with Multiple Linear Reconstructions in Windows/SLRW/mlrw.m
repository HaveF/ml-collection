
%    This code is written by Shiming Xiang
% 
%    Address:    National Laboratory of Pattern Recognition (NLPR), 
%                      Institute of Automation, Chinese Academy of Sciences, Beijing, 100190, China
%    Email:        smxiang@gmail.com

%    This software can be used freely for research purposes.
%    Published reports of research  using this code (or a modified version) should cite the 
%    article that describes the algorithm: 
%
%
%    Shiming Xiang, Chunhong Pan, Feiping Nie, Changshui Zhang: 
%    Interactive Image Segmentation With Multiple Linear Reconstructions in Windows. 
%    IEEE Transactions on Multimedia, Vol. 13, No.2, pp. 342-352, 2011.

%    This code has been compiled and tested using matlab 6.5  and matlab     7.0

%    version 1.0, 
%    updated by Shiming Xiang,  Oct. 19, 2011  

%==========================================================================

%==========================================================================
% ABOUT COMPUTATION TIME:
%
% With pure matlab codes, finsihing all the computaions may cost hundreds  of seconds. 
% For  example,  Finishing all of the pixels related to senetence  "M(ind, ind) = M(ind, ind) + ipl" in function
% "mlrw_calcu_global_matrix.m" may cost hundreds of seconds on PC. However, this step of adding local matrix
% to global matrix can be easily impleemnted with C language (For example, it  may cost less than 1 sencod 
% for images with 400 * 400 pixels on most PCs.).  
%
% We also wrote some C codes to replace some of the-time-consuming steps. But, for simplicity and 
% easy readability, here we only include the pure matlab codes of the algorithm.
%==========================================================================


function Y = mlrw(Img, Img_labeled,  mask_color, options)

% Img:                            the source image to be segmented, RGB color
% Img_labeled:               the labeled image to be segmented, RGB color image, uncompressed 
% mask_color:                used to extract the user specified foreground and background pixels
% options:
%         --- type_interaction:                tyep of interaction
%                                                       = 1;  interactive image segmentation (return hard edge)
%                                                       = 0;  interactive image matting (return soft edge)
%         --- gamma_global:                  the global regularization parameter,  gamma in Eq.(17)
%         --- lambda_local:                    the local regularization  parameter, lambda in Eq.(9)

% --------------------------------------------------------------------------------------------------------
% parameter check
% --------------------------------------------------------------------------------------------------------
if nargin <  4
   options.type_interaction = 1; 
   options.gamma_global     = 10000.0;
   options.lambda_local       = 0.0001;   
end

% check again as some fields can  also be missed even in the case that "options" is transferred
type_interaction = 1;
if isfield(options, 'type_interaction')
    type_interaction = options.type_interaction;
end

gamma_global = 10000.0;
if isfield(options, 'gamma_global')
    gamma_global = options.gamma_global;
end

lambda_local = 0.0001;
if isfield(options, 'lambda_local')
    lambda_local = options.lambda_local;
end

% --------------------------------------------------------------------------------------------------------
% convert image to matrix
% --------------------------------------------------------------------------------------------------------
X = cvrt_img_to_matrix(Img);                    %scan the image row by row. In X, each column is a data point, which is normalized by 255 (divided by 255) 
[RR, CC] = size(X);
X = X + 0.00000001 * randn(RR, CC);      % this sentence is actually unnecessary although I maintain it here

% --------------------------------------------------------------------------------------------------------
% extract the user labeled pixels
% --------------------------------------------------------------------------------------------------------
[height, width] = size(Img(:,:,1));
N = height * width;                                                         % number of pixels
Inds = extract_label_from_image(Img_labeled, mask_color);  % record them, row-by-row;

% --------------------------------------------------------------------------------------------------------
% constrcu the inforamtion related to the user labeled pixels
% --------------------------------------------------------------------------------------------------------
label_foreground = Inds{1, 1};
label_background = Inds{2, 1};
index_labeled = [label_foreground; label_background];
Y0 = zeros(N,1);

if type_interaction == 1                                         % interactive image segmentation
    Y0(label_foreground) = 1;
    Y0(label_background) = -1;
else
    Y0(label_foreground) = 1;                                % interactive image matting
    Y0(label_background) = 0;                               % note that although  we do not report the matting results  in the paper, it can also do matting   
end

clear Inds;
clear label_foreground;
clear label_background;

% --------------------------------------------------------------------------------------------------------
% construct the spatial neighborhood
% --------------------------------------------------------------------------------------------------------
win_width = 3;     % note that it is unnecessary to give a larger number. 3 is enough to obtain good results 
win_height = 3;
nb = calcu_only_spatial_neighorhood(width, height, win_width, win_height);     % including itself

% --------------------------------------------------------------------------------------------------------
% construct the global error matrix 
% --------------------------------------------------------------------------------------------------------
 
M = mlrw_calcu_global_matrix(X, nb, lambda_local);    % very slow!!!!!!!, commented by Xiang, Oct. 19, 2011
                                                                                     % Actually, this function can be totally replaced by C language


clear X;
clear nb;

% --------------------------------------------------------------------------------------------------------
% solve the linear system
% --------------------------------------------------------------------------------------------------------
Y0 = gamma_global .* Y0;
M = mlrw_regularized_matrix(M, index_labeled, gamma_global);
M = sparse(M);
Y = M \ Y0;
Y = reshape(Y, width, height)';

if type_interaction == 1                                         % interactive image segmentation
    ind = find(Y <= 0); 
    Y(ind) = 0; 
    ind = find(Y > 0); 
    Y(ind) = 255;
end

return;

    




