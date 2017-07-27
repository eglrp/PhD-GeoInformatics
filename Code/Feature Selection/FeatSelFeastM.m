%FeatSelFeastM Trainable mapping for FEAST library feature selection algorithms
%
%   [W,R] = FeatSelFeastM(A, CRIT, K, feastCell)
%   [W,R] = A*FeatSelFeastM([], CRIT, K, feastCell)
%   [W,R] = A*FeatSelFeastM(CRIT, K, feastCell)
%   [W,R] = FeatSelFeastM(A, CRIT, K, feastCell)
%   [W,R] = A*FeatSelFeastM([], CRIT, K, feastCell)
%   [W,R] = A*FeatSelFeastM(CRIT, K, feastCell)
%
% INPUT
%   A    Training dataset
%   CRIT Name of the FEAST algoritm see "help feast"
%   K    Number of features to select (default: sort all features)
%   feastCell  Cell array of additional args for FEAST (optional)
%
% OUTPUT
%   W    Feature selection mapping
%   R    Matrix with criterion values
%
% DESCRIPTION
% Wrapper for FEAST routines - see "help feast" for options.  
% For K = 0 all features are
% selected, but reordered according to the criterion. The result W can be
% used for selecting features using B*W.
% The selected features are stored in W.DATA and can be found by +W.
% In R, the search is reported step by step as:
%
% 	R(:,1) : number of features
% 	R(:,2) : criterion value
% 	R(:,3) : added / deleted feature
%
function [w, r] = FeatSelFeastM(varargin)

varargin = shiftargin(varargin, {'char', 'double'});
argin = setdefaults(varargin, [], 'jmi', 0, {});
if mapping_task(argin, 'definition')
    w = define_mapping(argin, 'untrained', 'FEAST feature selection wrapper');
    return
end

[a, crit, ksel, feastCell] = deal(argin{:});

[m, k, c] = getsize(a);
featlist = getfeatlab(a);

% If KSEL is not given, return all features.
if (ksel == 0)
    ksel = k;
end

isvaldfile(a, 1, 2); % at least 1 object per class, 2 classes
a = testdatasize(a);

% if strcmpi(crit, 'jmi')
%     ranked = JMIC(a, ksel, feastCell{:});
% else
ranked = feast(crit, ksel, +a, getnlab(a), feastCell{:});
% end

% There is no criterion value ???

r = [[1:ksel]', zeros(ksel, 1), ranked(:)];
J = ranked';
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
w = setname(w, ['FEAST ' crit ' feature selection']);

return
