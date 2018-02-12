%FEATSELI Trainable mapping for feature clustering and ranking
%
%   [W,R] = FeatSelLassoM(A, LAMBDA, K, T)
%   [W,R] = A*FeatSelLassoM([], LAMBDA, K, T)
%   [W,R] = A*FeatSelLassoM(LAMBDA, K, T)
%   [W,R] = FeatSelLassoM(A, LAMBDA, K, N)
%   [W,R] = A*FeatSelLassoM([], LAMBDA, K, N)
%   [W,R] = A*FeatSelLassoM(LAMBDA, K, N)
%
% INPUT
%   A    Training dataset
%   LAMBDA Regularisation parameter (larher = less features)
%        (default: 100)
%   K    Number of features to select (default: sort all features)
%   T    Tuning dataset (optional)
%   N    Number of cross-validation folds (optional)
%
% OUTPUT
%   W    Feature selection mapping
%   R    Matrix with criterion values
%
% DESCRIPTION
% LASSO feature selection. LAMBDA sets the regularisation parameter.
% criterion used by the feature evaluation routine FEATEVAL. If the dataset
% T is given, it is used as test set for FEATEVAL. For K = 0 all features are
% selected, but reordered according to the weights. The result W can be
% used for selecting features using B*W.
% The selected features are stored in W.DATA and can be found by +W.
% In R, the search is reported step by step as:
%
% 	R(:,1) : number of features
% 	R(:,2) : criterion value
% 	R(:,3) : added / deleted feature
%
function [w, r] = FeatSelLassoM(varargin)

varargin = shiftargin(varargin, {'char', 'prmapping'});
argin = setdefaults(varargin, [], 100, 0, []);
if mapping_task(argin, 'definition')
    w = define_mapping(argin, 'untrained', 'LASSO feature selection');
    return
end
fcrCell = {};
if length(argin) > 4
    fcrCell = argin(5:end);
    argin(5:end) = [];
end
[a, lambda, ksel, t] = deal(argin{:});

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

% convert the data set to target data set (assuming it is crisp label)
if ~strcmpi(getlabtype(a), 'targets')
    nl = getnlab(a);
    tgta = prdataset(+a);
    tgta = setlabtype(tgta, 'targets');
    tgta = settargets(tgta, nl);
else
    tgta = a;
end

tgta = tgta * scalem(tgta, 'variance'); % crucial for weights to be comparable

wlasso = lassor(tgta, lambda);

% NB the first coeff in w is the constant/offset and must be discarded
weights = +wlasso;
weights = weights((2:end));

fprintf('LASSO selected %d features of %d.  %d requested\n', sum(weights > 1e-9), k, ksel)
if ksel < k && sum(weights > 1e-9) < ksel
    error('Not enough features were selected - decrease k or increase lambda')
end

[weightSort, featIdx] = sort(-abs(weights(:)));
r = [[1:k]', abs(weightSort(:)), featIdx(:)];

featIdx = featIdx(weights > 1e-9);
J = featIdx(1:min(length(featIdx), ksel))';

% Return the mapping found.
w = featsel(k, J);
w = setmapping_type(w, 'trained');
w = setsize(w, [k, length(J)]);
if ~isempty(featlist)
    w = setlabels(w, featlist(J, :));
end
w = setname(w, 'LASSO FeatSel');

return
