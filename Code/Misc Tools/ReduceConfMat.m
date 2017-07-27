function [CmR] = ReduceConfMat(Cm, TgtIdx, isNorm)
%TgtIdx is cell array of indices for each class
if (isempty(TgtIdx))
    TgtIdx = mat2cell(1:size(Cm, 2), 1, ones(size(Cm, 2),1)); %no reduction
end
if (nargin < 3)
    isNorm = max(Cm(:)) <= 1;
end
numc = size(Cm, 1);
allIdx = 1:numc;

if ~iscell(TgtIdx) %old 2 class version
    bgIdx = setdiff(allIdx, TgtIdx);

    CmR(1, 1) = sum(sum(Cm(bgIdx, bgIdx), 2), 1);
    CmR(1, 2) = sum(sum(Cm(bgIdx, TgtIdx), 2), 1);
    CmR(2, 1) = sum(sum(Cm(TgtIdx, bgIdx), 2), 1);
    CmR(2, 2) = sum(sum(Cm(TgtIdx, TgtIdx), 2), 1);

    if (isNorm)
        n = [length(bgIdx), length(bgIdx); length(TgtIdx) length(TgtIdx)]; %wtf? - this looks rubbish.  Actually it is to weight each class equally not according to priors.
        %I THINK THIS SHOULD BE n' - CHECK
        CmR = CmR./n; %this will only work if the num of objects in each class is == or if you want to treat the data in this way i.e. as if it was equal priors
    end
else
    if (length(TgtIdx) == 1) %2 classes, infer the other one
        TgtIdx{2} = setdiff(allIdx, TgtIdx);
    end

    for i = 1:length(TgtIdx)
        for j = 1:length(TgtIdx)
            CmR(i, j) = sum(sum(Cm(TgtIdx{i}, TgtIdx{j}), 2), 1);
        end
        n(1,i) = length(TgtIdx{i});
    end

    if (isNorm) %this is the part where 'prior independance' is achieved
%         warning('Weighting as equal priors'); 
        CmR = CmR./repmat(n(:), 1, length(n)); %NB this will work as if the num of objects in each unlumped class is ==
    end
end
