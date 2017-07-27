function [W] = opencvdtreec(a, depth, optionCell)

% 	checktoolbox('libsvm');
if nargin < 3
    optionCell = {};
end
if nargin < 2
    depth = 10;
end

if nargin < 1 | isempty(a)
    
    W = prmapping(mfilename, {depth, optionCell});
    W = setname(W, 'opencvdtreec');
    return;
elseif ~ismapping(depth)
    islabtype(a,'crisp','soft');
    isvaldfile(a,1,2); % at least 1 object per class, 2 classes
    a = testdatasize(a,'objects');
    [m,k,c] = getsize(a);
    nlab = getnlab(a);
    
    if isempty(optionCell)
        optionCell = {'Priors', getprior(a), 'MaxDepth', depth, 'Use1seRule', false, ...
        'UseSurrogates', false, 'CVFolds', 5, 'TruncatePrunedTree', true, 'MinSampleCount', min(classsizes(a))/50};
    end
    
    clfr = cv.DTree;
    clfr.train(+a, nlab, optionCell{:});
    lablist = getlablist(a);
    
    W = prmapping(mfilename,'trained', {clfr}, lablist, size(a,2), c);
    
    W = setname(W, 'opencvdtreec');
    W = setcost(W, a);
else
    v = depth;
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
