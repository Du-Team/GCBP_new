function  C = mergeCoreNd(distLast,coreNd)
% Clustering the core nodes based on the distance relationship between them

    % Find out the proximity of the distances between the core points
    distCore = distLast(coreNd,coreNd);
    link = (distCore==1);
    C = adj2cluster(link);
    a=1;
    for i=1:length(C)
        C{i} =coreNd(C{i})';
    end

end