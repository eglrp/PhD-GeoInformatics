function res = FsStabilityEval(res, varargin)
%FsStabilityEval Evaluates feature selection stability on output of BootstrapFsEval 
%   Evaluates feature stability

ModifyDefaultArgs(varargin)
numBootStraps = length(res.FeatIdx);
% find the Pearson correlation coefficient between importance scores if
% they exist.  Note this only really works for ranking / featseli approaches and
% not for fs, be etc.  Finds on full set of features.
res.PearsonCC = corrcoef(res.FeatScore);

% find Spearman's rank correlation coefficient if there are importances for
% sorting.  Note one needs the proper ranks for all features i.e. you will
% need to add/remove all terms for fs/be respectively.  Finds on full set of features.
res.SpearmanRCC = corr(res.FeatRank, 'type', 'Spearman');

% find Tanimoto distance between all pairs of feature index sets. Finds on
% best numFeatures features
if ~iscell(res.FeatIdx)
    %res.FeatIdx = mat2cell(res.FeatIdx, size(res.FeatIdx, 1), ones(1, size(res.FeatIdx, 2)));
    res.FeatIdx = mat2cell(res.FeatIdx', ones(1, size(res.FeatIdx, 2)), size(res.FeatIdx, 1));
end
res.Tanimoto = zeros(numBootStraps, numBootStraps);
for i = 1:numBootStraps
    for j = 1:numBootStraps
        %D(:, j) = sum(abs(A - repmat(+B(j, :), m, 1)), 2);
        res.Tanimoto(i, j) = 1 - ((length(res.FeatIdx{i}) + length(res.FeatIdx{j}) - 2 * length(intersect(res.FeatIdx{i}, res.FeatIdx{j})))/...
            (length(res.FeatIdx{i}) + length(res.FeatIdx{j}) - length(intersect(res.FeatIdx{i}, res.FeatIdx{j}))));
    end
end

ll = cellfun(@length, res.FeatIdx);

% find consistency (see Brown et al 2012)
Dc = zeros(numBootStraps, numBootStraps);
if all(ll==ll(1)) %numFeatures > 0
    for i = 1:numBootStraps
        for j = 1:numBootStraps
            %D(:, j) = sum(abs(A - repmat(+B(j, :), m, 1)), 2);
            %this assumes length(res.FeatIdx{i}) == length(res.FeatIdx{j}) which may
            %not always be the case but will be with FEAST and where
            %numFeatures is spec'd
            r = length(intersect(res.FeatIdx{i}, res.FeatIdx{j}));
            k = length(res.FeatIdx{i});
            n = size(res.FeatScore, 1);  % check?
            Dc(i, j) = (r*n - k*k)/(k*(n - k));
        end
    end
end

if all(ll==ll(1))
    res.FeatIdx = cell2mat(res.FeatIdx')';
else
    res.FeatIdx = res.FeatIdx;
end

pDs = triu(res.Tanimoto, 1);
res.TanimotoStability = sum(pDs(:))/sum(1:numBootStraps-1);
pDw = triu(res.PearsonCC, 1);
res.PearsonCorrCoeffStab = sum(pDw(:))/sum(1:numBootStraps-1);
pDr = triu(res.SpearmanRCC, 1);
res.SpearmanRankCorrCoeffStab = sum(pDr(:))/sum(1:numBootStraps-1);
pDc = triu(Dc, 1);
res.Consitency = sum(pDc(:))/sum(1:numBootStraps-1);

