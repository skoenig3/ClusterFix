function [silh] = InterVSIntraDist(X, clust)
%Iter vs Intra distance points by cluster-mod of SILHOUETTE (just doesn't plot
%by Seth Koenig October 21, 2012
[idx,cnames] = grp2idx(clust);
n = length(idx);
k = length(cnames);
count = histc(idx(:)',1:k);
mbrs = (repmat(1:k,n,1) == repmat(idx,1,k));

% Get avg distance from every point to all (other) points in each cluster
myinf = zeros(1,1,class(X));
myinf(1) = Inf;
avgDWithin = repmat(myinf, n, 1);
avgDBetween = repmat(myinf, n, k);
for j = 1:n
    distj = sum((X - X(repmat(j,n,1),:)).^2, 2);
    % Compute average distance by cluster number
    for i = 1:k
        if i == idx(j)
            avgDWithin(j) = sum(distj(mbrs(:,i))) ./ max(count(i)-1, 1);
        else
            avgDBetween(j,i) = sum(distj(mbrs(:,i))) ./ count(i);
        end
    end
end

% Calculate the silhouette values
minavgDBetween = min(avgDBetween, [], 2);
silh = (minavgDBetween - avgDWithin) ./ max(avgDWithin,minavgDBetween);
s = silh;
end