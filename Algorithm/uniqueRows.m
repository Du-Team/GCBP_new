function [c,indA,indC] = uniqueRows(a) 
% Real nodes used for computation

    numRows = size(a,1);
    [sortA,indSortA] = sortrows(a);
    groupsSortA = sortA(1:numRows-1,:) ~= sortA(2:numRows,:);
    groupsSortA = any(groupsSortA,2);
    groupsSortA = [true; groupsSortA];
    groupsSortA = full(groupsSortA);
    c = sortA(groupsSortA,:); 
    indA = indSortA(groupsSortA);
    indCTemp = cumsum(groupsSortA);  
    indC(indSortA) = indCTemp;  
    indC = indC';
end
