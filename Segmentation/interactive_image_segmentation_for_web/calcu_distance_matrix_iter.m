%   version 2.0 --Jun/2010 
%   version 1.0 --Feb/2008 

%   written by Shiming Xiang ( smxiang@gmail.com ) 


function dist = calcu_distance_matrix_iter(X)

% X:     the input matrix, each column is a data point;

% dist:  squared distance: size N * N 

[D, N] = size(X);
X2 = sum(X.^2, 1);              % square each element, and then add each column

dist =  zeros(N, N);

for i = 1: N
    XX = X(:,i);                % extract the i-th column
        
    XX2 = (XX' * XX) * ones(1, N);
    dist_temp = XX2 - 2 * XX' * X;
    distance = X2 + dist_temp;
    
    dist(i, :) = distance;
end



