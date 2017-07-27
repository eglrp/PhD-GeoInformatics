function [W] = opencvrtreec(a, dum, optionCell)

% clfr.train(+tr, getnlab(tr), 'Priors', [1 1 1]/3, 'MaxNumOfTreesInTheForest', 10, 'NActiveVars', 4, 'CalcVarImportance', true, 'MaxDepth', 5, 'ForestAccuracy', 0.05);
if nargin < 3
    optionCell = {}; %'Priors', ones(1, getsize(a, 3))/getsize(a, 3)};
end
if nargin < 2
    dum = [];
end

if nargin < 1 | isempty(a)
    
    W = prmapping(mfilename, {dum, optionCell});
    W = setname(W, 'opencvrtreec');
    return;
elseif ~ismapping(dum)
    islabtype(a, 'crisp', 'soft');
    isvaldfile(a, 1, 2); % at least 1 object per class, 2 classes
    a = testdatasize(a, 'objects');
    [m, k, c] = getsize(a);
    nlab = getnlab(a);

    if isempty(optionCell)
        optionCell = {'Priors', ones(1, getsize(a, 3))/getsize(a, 3)};
    end

    clfr = cv.RTrees;
    clfr.train(+a, nlab, optionCell{:});
    lablist = getlablist(a);

    W = prmapping(mfilename, 'trained', {clfr}, lablist, size(a,2), c);

    W = setname(W, 'opencvrtreec');
    W = setcost(W, a);
else
    v = dum;
    clfr = v.data{1};
    response = clfr.predict(+a);
    out = zeros(size(a, 1), size(v, 2));
    
    %         out = zeros(length(yhat),getsize(a,3));
    for c = 1:size(v, 2)
        out(response==c, c) = 1;
    end
    %         out = vote./model.ntree;
    W = setdat(a, out, v);
end

return;
