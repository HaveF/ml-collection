function [label] = mykmeans(X,k)
    figure(1);
    plot(X(:,1),X(:,2),'o');
    colors = 'bgrcmyk';
    clustercolors = 'rcmykbg';
    hold on;
    [s1 s2] = size(X);
    xmin = min(X);
    xmax = max(X);
    label = zeros(1,s1);
    clusters = X(randsample(s1,k),:);

    hc= plot(clusters(:,1),clusters(:,2),'.','MarkerSize', 40);
    
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
        
        plot(X(idx,1),X(idx,2),[colors(i) 'o']);
        
        
    end
    delete(findobj(hc,'Type','Line'));
    hc = plot(clusters(:,1),clusters(:,2),'k.','MarkerSize', 40);
    delta = norm(clusters - oldclusters)
    pause(1)
    end
    
end