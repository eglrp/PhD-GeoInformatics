%FEATSELI Trainable mapping for feature clustering and ranking
%
%   [W,R] = FeatSelMultiViewM(A, LAMBDA, K, T)
%   [W,R] = A*FeatSelMultiViewM([], LAMBDA, K, T)
%   [W,R] = A*FeatSelMultiViewM(LAMBDA, K, T)
%   [W,R] = FeatSelMultiViewM(A, LAMBDA, K, N)
%   [W,R] = A*FeatSelMultiViewM([], LAMBDA, K, N)
%   [W,R] = A*FeatSelMultiViewM(LAMBDA, K, N)
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
% Multiview feature selection as in Chen et al 2017
% 	R(:,1) : number of features
% 	R(:,2) : criterion value
% 	R(:,3) : added / deleted feature
%
function [w, r] = FeatSelMultiViewM(varargin)

varargin = shiftargin(varargin, {'char', 'prmapping'});
argin = setdefaults(varargin, [], naivebc, 0, []);
if mapping_task(argin, 'definition')
    w = define_mapping(argin, 'untrained', 'Multiview Featsel');
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

fprintf('Affinity propagation:\n')
% Note that the scaling below makes a huge difference to clustering
% performance.  Scaling to unit variance works much better, but the
% paper only talks of scaling to range [0,1].  Their dist similarity measure
% is a weakness as closely spaced features are not the only redundant
% ones, features that are linearly scaled versions of each will be
% redundant too but they are not close.  Then there are non-lin dependencies too.
%     dataNorm = (data*scalem(data, 'domain'));  % according to chen et al 2017
dataNorm = (a * scalem(a, 'variance'));
S = -(+dataNorm)' * proxm((+dataNorm)', 'distance', 2);
%     S = -distm((+dataNorm)'); % -ve euclidean distance betw feats
m = size(S, 1); % num feats
tmp = triu(S, 1) + tril(S, -1); % kind of unnecessary as Sii = 0 already
pref = 1.1 * sum(tmp(:)) / (m * (m - 1)); % from paper but inc slightly otherwise we dont get enough features

[idx, netsim, i, unconverged, dpsim, expref] = apcluster(S, pref);
if unconverged
    S = S + 1e-9 * randn(size(S, 1), size(S, 2));
    [idx, netsim, i, unconverged, dpsim, expref] = apcluster(S, pref, 'maxits', 1000, 'dampfact', 0.8);
    if unconverged
        error('ERROR: apcluster unconverged');
    end
end

nclust = length(unique(idx));
tmp(unique(idx)) = 1:nclust;
clustlab = tmp(idx); % labels 1 indexed
exemplars = unique(idx);
fprintf('Number of clusters: %d\n', nclust);

count = 1;
for i = 1:nclust
    clustFeatIdx{i} = find(clustlab == i);
    clustFeatCount(i) = length(clustFeatIdx{i});
    fl = cellstr(getfeatlab(a));
    fprintf('Cluster %d, Exemplar %sf\n', i, fl{exemplars(i)});
    fprintf('%s, ', fl{clustFeatIdx{i}});
    fprintf('\n');
end

%% set up the matrices for quadprog as per the paper and my comments in Mendeley
[n, m] = size(a); %beware of redef
if false  % excl ||E||
    H = spalloc(2*(n + m), 2*(n + m), sum(clustFeatCount.^2));  %zeros(2*(n + m), 2*(n + m));
    clustOne = spalloc(n + m, n + m, sum(clustFeatCount.^2));
else  % incl ||E||
    H = spalloc(2*(n + m), 2*(n + m), sum(clustFeatCount.^2) + n);  %zeros(2*(n + m), 2*(n + m));
    clustOne = spalloc(n + m, n + m, sum(clustFeatCount.^2) + n);
    clustOne(m+1:end, m+1:end) = speye(n); %I think the paper has an issue here - 
                                           %this should give an L2 norm for
                                           %E, I can't figure out a way to
                                           %get an L1,2 norm as the group
                                           %components of e are hidden
    %clustOne(:, m+1:end) = 1; %I'm not sure of this but it is for E
end
%H = spalloc(2*(n + m), 2*(n + m), (n + m)*(n + m));  %zeros(2*(n + m), 2*(n + m));
for i = 1:nclust
    for j = 1:length(clustFeatIdx{i})
        clustOne(clustFeatIdx{i}(j), clustFeatIdx{i}) = 1;
    end
end

H(n+m+1:end, n+m+1:end) = clustOne;

% (n x 2(m+n)) . 2(m+n)
Aeq = spalloc(n, 2*(m+n), m*n + n);
% domain scaling is what they spec in the paper although it seems a poor
% choice.  outliers would do weird things.
Aeq = [+(a * scalem(a, 'domain')), lambda * speye(n), spalloc(n, m+n, 0)];
%Aeq = [+a, lambda * speye(n), spalloc(n, m+n, 0)];
beq = [getnlab(a)];
beq = beq - mean(beq);

A = [speye(m+n), -speye(m+n); -speye(m+n), -speye(m+n)];
b = zeros(2*(m + n), 1);

lb = zeros(2*(m + n), 1);
lb(1:m+n) = -inf;
lb(m+n+1:end) = 0;

options = optimoptions('quadprog');
% options.Algorithm = 'trust-region-reflective';

[x, fval, exitflag, output, lambda] = quadprog(H, [], A, b, Aeq, beq, lb, [], [], options);

[xsort, sortIdx] = sort(-abs(x(1:m)));
% x(sortIdx(1:ksel))
% fl(sortIdx(1:ksel))

% Sort the features by criterion value (maximum first).
% note that last column is rank and not idx as in other feat selectors
r = [[1:k]', x(sortIdx(1:k)), sortIdx(1:k)];
J = sortIdx(1:ksel)';
%
% [critval_sorted, J] = sort(-critval);
% 	r = [[1:k]', -critval_sorted, J];
% 	J = J(1:ksel)';


% Return the mapping found.
w = featsel(k, J);
w = setmapping_type(w, 'trained');
w = setsize(w, [k, length(J)]);
if ~isempty(featlist)
    w = setlabels(w, featlist(J, :));
end
w = setname(w, 'MultiView FeatSel');

return
