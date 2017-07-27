%prmapping wrapper for R backend random forest
function [W,model] = librandomforestc(a,ntree,mtry,extra_options)

% 	checktoolbox('libsvm');
    if nargin < 4
        extra_options = [];
    end
    if nargin < 3
        mtry = 0;
    end
    if nargin < 2
        ntree = 0;
    end
	
	if nargin < 1 | isempty(a)
		W = prmapping(mfilename,{ntree,mtry,extra_options});
		W = setname(W,'librandomtreec');
		return;
    elseif ~ismapping(ntree)
		islabtype(a,'crisp','soft');
		isvaldfile(a,1,2); % at least 1 object per class, 2 classes
		a = testdatasize(a,'objects');
		[m,k,c] = getsize(a);
		nlab = getnlab(a);
        extra_options.classwt = getprior(a);
        extra_options.importance = true; %assess feature importance
        extra_options.do_trace = 0;
        model = classRF_train(+a,nlab, ntree, mtry, extra_options);
		lablist = getlablist(a);
        
		W = prmapping(mfilename,'trained', {model}, lablist(model.new_labels, :), size(a,2), c);
		
		W = setname(W,'librandomtreec');
		W = setcost(W,a);
    else
        v = ntree;
        model = v.data{1};
        [yhat, vote] = classRF_predict(+a, model);
%         out = zeros(length(yhat),getsize(a,3));
%         for c = 1:getsize(a,3)
%             out(yhat==c, c) = 1;
%         end
        out = vote./model.ntree;
        W = setdat(a, out, v);
    end
	
return;
