%   version 2.0 --Jun/2010 
%   version 1.0 --Feb/2008 

%   written by Shiming Xiang ( smxiang@gmail.com ) 


function LL = construct_coefficient_matrix(X)
% X:  each column is a data point, 

%return
% LL:    return the coefficient matrix


%    all the data points are fully used to construct the spline

[dim, num] = size(X);          % the size of the data matrix      

distance = calcu_distance_matrix_iter(X);          % note that:     returned squarded distances

% adopt the spline
lambda = 0.0001;
distance   =  distance .* log(distance + 0.0001) + lambda * eye(num);    % there exists a scale of 2.0 with log function, this will not change the results as alpha can be viewed as a new alpha =  (2.0*alpha)

%construct the coefficient matrix
P = ones(num, 1);
P = [P, X'];         
C = zeros(dim + 1, dim + 1);
LL = [distance, P; P', C];

return