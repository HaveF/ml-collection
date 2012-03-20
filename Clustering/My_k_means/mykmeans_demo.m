function mykmeans_demo()
load('data.mat');
mykmeans(x',4);
end


function [label] = mykmeans(X,k)
    plot(X(:,1),X(:,2),'o');
    colors = 'bgrcmyk';
    hold on;
    [s1 s2] = size(X);
    xmin = min(X);
    xmax = max(X);
    label = zeros(1,s1);
    clusters = zeros(k,s2);
    for i= 1:k
        for j = 1:s2
            clusters(i,j) = xmin(j)+(xmax(j)-xmin(j))*rand();
        end
        hc(i) = plot(clusters(i,1),clusters(i,2),[colors(8-i) '.'],'MarkerSize', 40);
    end
    delta = 100;
    ditance = zeros(1,k);
    while delta>1e-8
    oldclusters = clusters;
    for i = 1:s1
        for j=1:k
            ditance(j) = norm(X(i,:)-oldclusters(j,:));
        end
        [iminvalue minidx ] = min(ditance);
        label(i) = minidx;
    end
    %update clusters
    for i=1:k
        idx = find(label==i);
        clusters(i,:) =  mean(X(idx,:));
        delete(findobj(hc(i),'Type','Line')); 
        plot(X(idx,1),X(idx,2),[colors(i) 'o']);
        hc(i) = plot(clusters(i,1),clusters(i,2),[colors(8-i) '.'],'MarkerSize', 40);
        
    end
    delta = norm(clusters - oldclusters)
    pause(1)
    end
    
end