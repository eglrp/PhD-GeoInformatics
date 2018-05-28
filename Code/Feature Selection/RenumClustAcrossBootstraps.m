function res = RenumClustAcrossBootstraps(res)
% renumbers clusters across bootstraps, so that identical clusters have
% same cluster number and unique clusters have unique cluster numbers

numBootStraps = size(res.FeatClustIdx, 2);
newClusterNum = max(res.FeatClustIdx(:))+1;
res.FeatClustIdxOrig =  res.FeatClustIdx;

% initialise renumClustList
for bi = 1:numBootStraps  % iner loop is clusters to renumber
    featClustIdx = res.FeatClustIdx(:, bi);
    clusterList = unique(featClustIdx);            
    renumClustList{bi} = zeros(1, length(clusterList));
end
    
for bi = 1:numBootStraps-1  % outer loop is clusters to compare to
    masterFeatClustIdx = res.FeatClustIdx(:, bi);
    masterClustList = unique(masterFeatClustIdx, 'stable');
    for bj = bi+1:numBootStraps  % iner loop is clusters to renumber
        compFeatClustIdx = res.FeatClustIdx(:, bj);
        compClusterList = unique(compFeatClustIdx, 'stable');            
%         renumClustList = zeros(1, length(compClusterList));
        for ci = 1:length(masterClustList)  % loop through master clusters
            masterClustFeatIdx = find(masterFeatClustIdx == masterClustList(ci));
            for cj = 1:length(compClusterList)  % loop through compare clusters
                compClustFeatIdx = find(compFeatClustIdx == compClusterList(cj));
                sd = setxor(masterClustFeatIdx, compClustFeatIdx);  % don't use setdiff - it is not symmetrical
                if isempty(sd)  %clusters are identical, give them the same number
                    renumClustList{bj}(cj) = masterClustList(ci);
                    break;
                end
            end
        end
        % assign new cluster numbers to those clusters that are unique
        numNewClusters = sum(renumClustList{bj} == 0);
        newnumClustList = renumClustList{bj};
        newnumClustList(renumClustList{bj} == 0) = newClusterNum:(newClusterNum + numNewClusters - 1);
        newClusterNum = newClusterNum + numNewClusters;

        renumFeatClustIdx = res.FeatClustIdx(:, bj);
        for cj = 1:length(compClusterList)
            idx = compFeatClustIdx == compClusterList(cj);
            renumFeatClustIdx(idx) = newnumClustList(cj);
        end
        res.FeatClustIdx(:, bj) = renumFeatClustIdx;
    end
end
% for bi = 2:numBootStraps  % iner loop is clusters to renumber
%     featClustIdx = res.FeatClustIdx(:, bi);
%     clusterList = unique(featClustIdx);            
%     renumFeatClustIdx = res.FeatClustIdx(:, bi);
%     for ci = 1:length(clusterList)
%         idx = featClustIdx == clusterList(ci);
%         renumFeatClustIdx(idx) = renumClustList{bi}(ci);
%     end
%     res.FeatClustIdx(:, bi) = renumFeatClustIdx;
% end
return
%% Orig 2016 Oct Code from BootstrapFsEval
    newClusterNum = length(res.ClustFeatNLab{1})+1;
    for b = 2:numBootStraps
        clustReNum = zeros(1, length(res.ClustFeatNLab{b}));
        for c1 = 1:length(res.ClustFeatNLab{1})
            for c = 1:length(res.ClustFeatNLab{b})
                sd = setdiff(res.ClustFeatNLab{1}{c1}, res.ClustFeatNLab{b}{c});
                if isempty(sd)  %clusters are identical, give them the same number
                    clustReNum(c) = c1;
                    break;
                end
            end
        end
        numNewClusters = sum(clustReNum == 0);
        clustReNum(clustReNum==0) = newClusterNum:(newClusterNum + numNewClusters - 1);
        newClusterNum = newClusterNum + numNewClusters;
        %this assumes that FeatClustIdx is numbered from 1 and
        %corresponding to order of ClustFeatNLab
        res.FeatClustIdx(:, b) = clustReNum(res.FeatClustIdx(:, b));
    end

