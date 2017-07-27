function [R] = NonDominantRank(X)
%As in Brown et al 2012 section 5.5
%   Using Pareto fronts, find rank of variables with > 1 fitness measure
%   Cols of x are the fitness measures, rows are objects / variables 
%   (eg cols are accuracy and stability, rows are method)
%   R is the rank of the cols of X

X_ = X;
R = zeros(size(X, 1), 1);
idx = 1:size(X, 1);
rank = 0;
for i = 1:size(X, 1)
    front = paretofront(X_);
    X_(front, :) = [];  %delete the dominant points
    R(idx(front)) = i;
    idx(front) = [];
    rank = rank + sum(front);
    if isempty(X_)
        break;
    end
end

end
