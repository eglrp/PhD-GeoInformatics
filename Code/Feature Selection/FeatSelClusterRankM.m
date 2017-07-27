%FEATSELI Trainable mapping for feature clustering and ranking
%
%   [W,R] = FeatSelClusterRankM(A, CRIT, K, T)
%   [W,R] = A*FeatSelClusterRankM([], CRIT, K, T)
%   [W,R] = A*FeatSelClusterRankM(CRIT, K, T)
%   [W,R] = FeatSelClusterRankM(A, CRIT, K, N)
%   [W,R] = A*FeatSelClusterRankM([], CRIT, K, N)
%   [W,R] = A*FeatSelClusterRankM(CRIT, K, N)
%
% INPUT
%   A    Training dataset
%   CRIT Name of the criterion or untrained mapping
%        (default: 'NN', i.e. the LOO 1-Nearest Neighbor error)
%   K    Number of features to select (default: sort all features)
%   T    Tuning dataset (optional)
%   N    Number of cross-validation folds (optional)
%
% OUTPUT
%   W    Feature selection mapping
%   R    Matrix with criterion values
%
% DESCRIPTION
% Feature clustering and ranking of K features using the dataset A. CRIT sets the
% criterion used by the feature evaluation routine FEATEVAL. If the dataset
% T is given, it is used as test set for FEATEVAL. For K = 0 all features are
% selected, but reordered according to the criterion. The result W can be
% used for selecting features using B*W.
% The selected features are stored in W.DATA and can be found by +W.
% In R, the search is reported step by step as:
%
% 	R(:,1) : number of features
% 	R(:,2) : criterion value
% 	R(:,3) : added / deleted feature
%
function [w, r] = FeatSelClusterRankM(varargin)

varargin = shiftargin(varargin, {'char', 'prmapping'});
argin = setdefaults(varargin, [], naivebc, 0, []);
if mapping_task(argin, 'definition')
    w = define_mapping(argin, 'untrained', 'Feature Clustering and Ranking');
    return
end
fcrCell = {};
if length(argin) > 4
    fcrCell = argin(5:end);
    argin(5:end) = [];
end
[a, crit, ksel, t] = deal(argin{:});

[m, k, c] = getsize(a);
featlist = getfeatlab(a);

% If KSEL is not given, return all features.
if (ksel == 0)
    ksel = k;
end

isvaldfile(a, 1, 2); % at least 1 object per class, 2 classes
a = testdatasize(a);
if isdataset(t)
    iscomdset(a, t);
end

% TO DO: pass criterion and other varargin
res = FeatureClusterRank(a, 'criterion', crit, fcrCell{:});
% select 1 feat from each cluster now
% [clustAcc_ clustIdx] = sort(res.ClustAcc);
% 
% fprintf('Cluster Ranking:\n----------------\n');
% feats = [];
% 
% % choose the most accurate feature from each cluster
% for i = 1:length(res.ClustAcc)
%     clusterFeatIdx = res.ClustFeatNLab{clustIdx(i)};
%     [featAcc_ sortClusterFeatIdx] = sort(res.FeatIndividualAcc(clusterFeatIdx));
%     feats(i) = clusterFeatIdx(sortClusterFeatIdx(1));
% end
% TO DO: allow user to select best feats

% Sort the features by criterion value (maximum first).
% note that last column is rank and not idx as in other feat selectors
r = [[1:k]', res.FeatClustScore(:), res.FeatClustRank(:)];
J = res.Feats(1:min(length(res.Feats), ksel))';
%
% [critval_sorted, J] = sort(-critval);
% 	r = [[1:k]', -critval_sorted, J];
% 	J = J(1:ksel)';


% Return the mapping found.
w = featsel(k, J);
w = setmapping_type(w, 'trained');
w = setsize(w, [k length(J)]);
if ~isempty(featlist)
    w = setlabels(w, featlist(J, :));
end
w = setname(w, 'Feature Clustering and Ranking');
w = setuser(w, res);

return
