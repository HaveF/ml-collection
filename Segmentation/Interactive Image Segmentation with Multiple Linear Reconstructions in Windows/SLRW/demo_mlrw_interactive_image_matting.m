
%    This code is written by Shiming Xiang
% 
%    Address:    National Laboratory of Pattern Recognition (NLPR), 
%                      Institute of Automation, Chinese Academy of Sciences, Beijing, 100190, China
%    Email:        smxiang@gmail.com

%    This software can be used freely for research purposes.
%    Published reports of research  using this code (or a modified version) should cite the 
%    article that describes the algorithm: 
%
%    Shiming Xiang, Chunhong Pan, Feiping Nie, Changshui Zhang: 
%    Interactive Image Segmentation With Multiple Linear Reconstructions in Windows. 
%    IEEE Transactions on Multimedia, Vol. 13, No.2, pp. 342-352, 2011.

%    This code has been compiled and tested using matlab 6.5  and matlab     7.0

%    version 1.0, 
%    updated by Shiming Xiang,  Oct. 19, 2011  

%==========================================================================



%==========================================================================
% WHAT'S NEW: !!!!!
% Note that: Our algorithm can also be extended for image matting, although we did not report this in the paper
% In this matlab package, the function for image matting with our method is supported.
% Please run the codes in file: demo_mlrw_interactive_image_matting.m
%==========================================================================


%==========================================================================
% ABOUT COMPUTATION TIME:
%
% With pure matlab codes, finishing all the computations may cost hundreds  of seconds. 
% For  example,  Finishing all of the pixels related to sentence  "M(ind, ind) = M(ind, ind) + ipl" in  function
% "mlrw_calcu_global_matrix.m" may cost hundreds of seconds on PC. However, this step of adding local matrix
% to global matrix can be easily implemented with C language.  (For example, it  may cost less than 1 sencod 
% for images with 400 * 400 pixels on most PCs.).  
%
% We also wrote some C codes to replace some of the-time-consuming steps. But, for simplicity and 
% easy readability, here we only include the pure matlab codes of the algorithm.
%==========================================================================


close all

% the source image to be segmented
file_index = '002';

current_work_folder  = pwd;          % Hope this matlab command can  work correctly
file_save = [current_work_folder, '\', file_index, '_matting.bmp'];


Img = imread( [file_index, '.bmp'] );
Img_labeled = imread( [file_index, '_label.bmp'] ); 

mask_color = [0,     255, 0;...       % foreground 
                       255, 0,     0];         % background

                   
% for interactive image matting, we must set  "options.type_interaction = 0;"    % must have  "~= 1" 
options.type_interaction = 0; 
Y = mlrw(Img, Img_labeled,  mask_color, options);   % parameter 'options' now  must be transferred for interactive imatge matting


%--------------------------------------------------------------------------
% % call it in other ways (and so on)  
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% % call it in other ways (and so on)  
% % 1 
% options.type_interaction = 0;                                                         % for image segmentaion,   set options.type_interaction = 1;
% Y = mlrw_for_interactive_image(Img, Img_labeled,  mask_color, options);

%--------------------------------------------------------------------------
% % 2 
% options.type_interaction = 0; 
% options.gamma_global     = 10000.0;
% Y = mlrw_for_interactive_image(Img, Img_labeled,  mask_color, options);

%--------------------------------------------------------------------------
% % 3 
% options.type_interaction = 0; 
% options.gamma_global     = 10000.0;
% options.lambda_local       = 0.0001;   
% Y = mlrw_for_interactive_image(Img, Img_labeled,  mask_color, options);

%--------------------------------------------------------------------------
% % 4 ...
%--------------------------------------------------------------------------


figure(1);
imshow(Y);
imwrite(Y,  file_save);

return;
