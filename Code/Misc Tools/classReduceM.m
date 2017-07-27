function [b] = classReduceM (a, classIdxCell)
%classIdxCell should be a cell array of source class indices to be
%combined into reduced classes

%this mapping does a class 'lumping' i.e. orig classes are grouped to
%produce a new reduced set of classes. 
%WARNING: usually one would want to this to treat the orig classes as =
%priors, to avoid biasing by these orig class priors in the lumped classes.
%Using this mapping will not achieve this in terms of performance measures
%however, only in terms of the training of the clfr


	prtrace(mfilename);

	% Default: perform PCA.
    if (nargin < 2)
        error('You must supply both arguments');
    end

    mapname = 'Class reduction mapping';
    flb = 1:length(classIdxCell);
    numc = 0;
    for i = 1:length(classIdxCell)
        numc = numc + length(classIdxCell{i});
    end

	% Empty mapping: return straightaway.
	if (isempty(a))
		b = mapping(mfilename, 'fixed', {classIdxCell}, flb(:), numc, length(classIdxCell));
		b = setname(b, mapname);
		return
	end

	%nodatafile(a);
	if ~isdataset(a)
		if isa(a,'double')
			a = dataset(a,1);   % make sure we have a dataset
		else
			error('nodatafile','Data should be given in a dataset or as doubles')
		end
	end

	islabtype(a,'crisp','soft');
% 	isvaldfile(a,1);   % at least 1 object per class

% 	[m, k, c] = getsize(a);

%     b = a(:, 1); %get id and class labels
    nl = getnlab(a);
    nlb = zeros(size(nl));
%     pa = getprior(a);
    for i = 1:length(classIdxCell)
        b(:, i) = max(+a(:, classIdxCell{i}), [], 2);
        pb(i) = 0;
        for j = 1:length(classIdxCell{i})
            nlb(nl == classIdxCell{i}(j)) = i;
%             pb(i) = pb(i) + pa(j);
        end
    end

    b = dataset(b);
    b = setfeatlab(b, flb(:));
%     , nlb(:)
    b = setlablist(b, [1:length(classIdxCell)]');
    b = setnlab(b, nlb(:));
    id = getident(a);
    if (~isempty(id))
        b = setident(b, id);
    end
%     b = setprior(b, pb); %??
     b = setprior(b, 0); %VERY dubious
%     b = setfeatlab(b, featLabCell); %relabel feats
%     w = setident(getident(b));
return
