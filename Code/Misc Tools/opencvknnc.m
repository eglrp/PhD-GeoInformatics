function [W] = opencvknnc(a, K, optionCell)

% 	checktoolbox('libsvm');
if nargin < 3
    optionCell = {};
end
if nargin < 2
    K = 3;
end

if nargin < 1 | isempty(a)
    
    W = prmapping(mfilename, {K, optionCell});
    W = setname(W, 'opencvknnc');
    return;
elseif ~ismapping(K)
    islabtype(a,'crisp','soft');
    isvaldfile(a,1,2); % at least 1 object per class, 2 classes
    a = testdatasize(a,'objects');
    [m,k,c] = getsize(a);
    nlab = getnlab(a);
    
    clfr = cv.KNearest(+a, nlab);
    lablist = getlablist(a);
    
    W = prmapping(mfilename, 'trained', {clfr, K}, lablist, size(a,2), c);
    
    W = setname(W, 'opencvknnc');
    W = setcost(W, a);
else
    v = K;
    K = v.data{2};
    clfr = v.data{1};
    response = clfr.predict(+a, 'K', K);
    out = zeros(size(a, 1), size(v, 2));
    for c = 1:size(v, 2)
        out(response == c, c) = 1;
    end
    W = setdat(a, out, v);
end

return;
