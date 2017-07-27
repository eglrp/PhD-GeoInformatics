%FeatSelReliefF Trainable mapping for RELIEFF feature selection
% wraps the stats toolbox function
%
%   [W,R] = FeatSelReliefF(A, K, KSEL)
%   [W,R] = A*FeatSelReliefF([], K, KSEL)
%   [W,R] = A*FeatSelReliefF(K, KSEL)
%
% INPUT
%   A    Training dataset
%   K    RELIEFF param
%   KSEL Num features to choose (optional)
%
% OUTPUT
%   W    Feature selection mapping
%   R    Matrix with criterion values
%
% DESCRIPTION
% RELIEFF using the dataset A. Wraps the stats toolbox function. 
% The result W can be used for selecting features using B*W.
% The selected features are stored in W.DATA and can be found by +W.
% In R, the search is reported step by step as:
%
% 	R(:,1) : number of features
% 	R(:,2) : weight
% 	R(:,3) : rank / index
%
function [w, r] = FeatSelReliefF(varargin)

varargin = shiftargin(varargin, {'double', 'double'});
argin = setdefaults(varargin, [], 10, 0);
if mapping_task(argin, 'definition')
    w = define_mapping(argin, 'untrained', 'RELIEFF');
    return
end
fcrCell = {};
if length(argin) > 3
    fcrCell = argin(4:end);
    argin(4:end) = [];
end
[a, kr, ksel] = deal(argin{:});

[m, k, c] = getsize(a);
featlist = getfeatlab(a);

% If KSEL is not given, return all features.
if (ksel == 0)
    ksel = size(a, 2);
end

isvaldfile(a, 1, 2); % at least 1 object per class, 2 classes
a = testdatasize(a);

[ranked, weight] = relieff(+a, getnlab(a), kr, 'prior', 'uniform', ...
    'method', 'classification', 'categoricalx', 'off');

% Sort the features by criterion value (maximum first).
[weightSort, featIdx] = sort(-weight(:));
r = [[1:k]', -weightSort, featIdx(:)];
J = ranked(1:min(k, ksel))';

% [tmp idx] = sort(-weight);
% assert(all(idx==ranked))

% Return the mapping found.
w = featsel(k, J);
w = setmapping_type(w, 'trained');
w = setsize(w, [k length(J)]);
if ~isempty(featlist)
    w = setlabels(w, featlist(J, :));
end
w = setname(w, 'RELIEFF');

return
