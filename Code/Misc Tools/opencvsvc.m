function [W] = opencvsvc(a, dum, optionCell)

% clfr.train((+tr), getnlab(tr), );
if nargin < 3
    %NB the below are defaults for my MSC
    optionCell = {'SVMType', 'C_SVC', 'KernelType', 'RBF', 'Gamma', 1, 'C', 1};
%     optionCell = {'SVMType', 'C_SVC', 'KernelType', 'RBF', 'Gamma', 5, 'C', 1, 'ClassWeights', double([1; 1; 1.5])};
end
if nargin < 2
    dum = [];
end

if nargin < 1 | isempty(a)
    
    W = prmapping(mfilename, {dum, optionCell});
    W = setname(W, 'opencvsvc');
    return;
elseif ~ismapping(dum)
    islabtype(a, 'crisp', 'soft');
    isvaldfile(a, 1, 2); % at least 1 object per class, 2 classes
    a = testdatasize(a, 'objects');
    [m, k, c] = getsize(a);
    nlab = getnlab(a);

    clfr = cv.SVM;
    clfr.train(+a, nlab, optionCell{:});
    lablist = getlablist(a);

    W = prmapping(mfilename, 'trained', {clfr}, lablist, size(a,2), c);

    W = setname(W, 'opencvsvc');
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
