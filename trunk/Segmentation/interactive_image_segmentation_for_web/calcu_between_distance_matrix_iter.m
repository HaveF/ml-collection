%   version 2.0 --Jun/2010 
%   version 1.0 --Feb/2008 

%   written by Shiming Xiang ( smxiang@gmail.com ) 


function dist = calcu_between_distance_matrix_iter(XR, XS)

% XR:     the reference data points, each column is a data point;
% XS:     the data points for calculating the distance matrix from XS to XR. Each column is a data point

% dist:  squared distance: each row is a distance vector for a data point w.r.t. XR




[D, Nr] = size(XR);
[D, Ns] = size(XS);

X2 = sum(XR.^2, 1);              % square each element, and then add each column
dist =  zeros(Ns, Nr);

% I do not thinkwritting the following sentences with a for-loop is a good selection. Actually we can try to use technique of block matrices to speed up the computaion.

for i = 1: Ns
    XX = XS(:,i);                                              % extract the i-th column
    XX2 = (XX' * XX) * ones(1, Nr);
    dist_temp = XX2 - 2 * XX' * XR;
    distance = X2 + dist_temp;
    dist(i, :) = distance;
end





