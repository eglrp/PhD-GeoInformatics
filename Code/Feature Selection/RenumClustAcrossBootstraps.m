function res = RenumClustAcrossBootstraps(res)
% renumbers clusters across bootstraps, so that identical clusters have
% same cluster number and unique clusters have unique cluster numbers

numBootStraps = size(res.FeatIdx, 2);
newClusterNum = max(res.FeatClustIdx(:))+1;
res.FeatClustIdxOrig =  res.FeatClustIdx;
for bi = 1:numBootStraps-1  % outer loop is clusters to compare to
    masterFeatClustIdx = res.FeatClustIdx(:, bi);
    masterClustList = unique(masterFeatClustIdx);
    for bj = bi+1:numBootStraps  % iner loop is clusters to renumber
        compFeatClustIdx = res.FeatClustIdx(:, bj);
        compClusterList = unique(compFeatClustIdx);            
        renumClustList = zeros(1, length(compClusterList));
        for ci = 1:length(masterClustList)  % loop through master clusters
            masterClustFeatIdx = find(masterFeatClustIdx == masterClustList(ci));
            for cj = 1:length(compClusterList)  % loop through compare clusters
                compClustFeatIdx = find(compFeatClustIdx == compClusterList(cj));
                sd = setdiff(masterClustFeatIdx, compClustFeatIdx);
                if isempty(sd)  %clusters are identical, give them the same number
                    renumClustList(cj) = masterClustList(ci);
                    break;
                end
            end
        end
        % assign new cluster numbers to those clusters that are unique
        numNewClusters = sum(renumClustList == 0);
        renumClustList(renumClustList == 0) = newClusterNum:(newClusterNum + numNewClusters - 1);
        newClusterNum = newClusterNum + numNewClusters;

        renumFeatClustIdx = res.FeatClustIdx(:, bj);
        for cj = 1:length(compClusterList)
            idx = compFeatClustIdx == compClusterList(cj);
            renumFeatClustIdx(idx) = renumClustList(cj);
        end
        res.FeatClustIdx(:, bj) = renumFeatClustIdx;
    end
end
