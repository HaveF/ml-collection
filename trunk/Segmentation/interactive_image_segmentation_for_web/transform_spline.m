%   version 2.0 --Jun/2010 
%   version 1.0 --Feb/2008 

%   written by Shiming Xiang ( smxiang@gmail.com ) 


function Y = transform_spline(CM, XR, XS)
%CM:    the coefficient matrix
%XR:     the reference points (used for training)
%XS:     the points to be transformed, each column is a data point (test data points)

%return
%Y:     the coordinate of the transformed data point, each column is a data point

%         transform Ys according to Yd. Ys are to be transformed: Y = f(Ys)  


[D, N] = size(XS);
dim = size(CM, 1);
Y = zeros(dim, N);

% actually, we can finish the latter calculations with a few sentences (see the commented sentences attached in the end of this file).
% But out-of-memory may occur for very large images. So we can divide the test data points into   small subsets. 
% Actually, in this way, the computation time will not be significantly increased

% here we divide the test data points for many subsets to avoid out-of-memory for very large images

num_pts_subgroup = 50000;    % if this number is also  large,  you can change it as 20000, 10000, or 5000
num_steps = floor(N / num_pts_subgroup);

for i = 1 : num_steps
    XS_temp = XS(:,  (i - 1) * num_pts_subgroup + 1 :  i * num_pts_subgroup);
    distance = calcu_between_distance_matrix_iter(XR, XS_temp);
    distance   =  distance .* log(distance + 0.0001);
    P = ones(1, num_pts_subgroup);
    P = [P; XS_temp];                                            % add the rows;
    P1 = [distance, P'];                
    Y(:, (i - 1) * num_pts_subgroup + 1 :  i * num_pts_subgroup) = CM * P1'; 
end

if  N - num_steps * num_pts_subgroup > 0
    XS_temp = XS(:,  num_steps * num_pts_subgroup + 1 :  N);
    distance = calcu_between_distance_matrix_iter(XR, XS_temp);
    distance   =  distance .* log(distance + 0.0001);
    num_pts_here = size(XS_temp, 2);
    P = ones(1, num_pts_here);
    P = [P; XS_temp];                                            % add the rows;
    P1 = [distance, P'];                
    Y(:, num_steps * num_pts_subgroup + 1 :  N) = CM * P1';     
end

return;




% %==============================================================
% % JUST FOR REFERENCE
% %==============================================================

% %==============================================================
% % the above codes can be simplified as follows (if the size of  image to be segmented  is small) 
% %==============================================================
% 
% [D, N] = size(XS);
% distance = calcu_between_distance_matrix_iter(XR, XS);    % out-of-memory problem may occur here 
% distance   =  distance .* log(distance + 0.0001);
% P = ones(1, N);
% P = [P; XS];                       % add the rows;
% 
% dim = size(CM, 1);
% Y = zeros(dim, N);
% P1 = [distance, P'];                                                              % out-of-memory problem may occur here          
% Y = CM * P1';
% return
%
% %==============================================================