%   version 2.0 --Jun/2010 
%   version 1.0 --Feb/2008 

%   written by Shiming Xiang ( smxiang@gmail.com ) 


function [X_all_mapped, CM] = spline_regression(X, Y, X_all) 
% X:                        each row is a data point;
% Y:                        the target values of X

% X_all:                    the test data points to be mapped

% return 
% X_all_mapped:    each column is a data point after mapped



CM = solve_linear_system(X', Y);              % Y = f(X)
 
% transform all data points;
X_all_mapped  = transform_spline(CM, X', X_all');        % for the test  data point,  each column is a data point after mapped
 

return;
