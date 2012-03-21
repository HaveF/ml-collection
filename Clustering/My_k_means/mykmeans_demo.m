function mykmeans_demo()
clc;
load('data.mat');
reslt = mygmm(x(:,1:100)',4);
%mykmeans(x',4);
end
