function [n2,dis]=uniform_sample_line(n,a)
dis=[0;cumsum(sqrt(sum((n(2:end,:)-n(1:end-1,:)).^2,2)))];
disn=linspace(dis(1),dis(end),a);
n2(:,1)= interp1(dis,n(:,1),disn);
n2(:,2)= interp1(dis,n(:,2),disn);
n2(:,3)= interp1(dis,n(:,3),disn);