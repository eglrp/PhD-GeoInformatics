function res = AffinityPropagation(data)
%Affinity propagation clustering as presented in X. Chen, G. Zhou, Y. Chen,
%G. Shao and Y. Gu (2017)
% data - prtools data set
data = data*scalem(data, 'variance');   % scale to unit variance
S = -distm((+data)');  % -ve euclidean distance betw feats
n = size(S, 1);     % num feats
tmp = triu(S, 1) + tril(S, -1);  % kind of unnecessary as Sii = 0 already
pref = sum(tmp(:))/(n*(n-1));   % from paper
S = tmp + diag(pref*ones(n, 1));

A = zeros(n, n);
R = zeros(n, n);
count = 0;
RAprev = R+A;

while true
    AS = A+S;
%     for r = 1:n
%         for c = 1:n
%             m1(c) = max(as([1:c-1; c+1:n]));
%         end
%         R(r, :) = S(r, :) - max(m1(idx));
%     end
    for c = 1:n
        R(:, c) = S(:, c) - max(AS(:, [1:c-1, c+1:n]), [], 2);
    end
    for r = 1:n
        Aineqj = min(0, diag(R) + sum(max(0, R([1:r-1, r+1:n], :)), 1)');
        Aieqj = sum(max(0, R([1:r-1, r+1:n], r)), 1);
        A(r, :) = Aineqj;
        A(r, r) = Aieqj;
    end
    if count > 100 || sum(sum(abs(RAprev - (R + A)))) < 1e-6
        break;
    end
    fprintf('.')
    RAprev = R+A;
    count = count + 1;
end
fprintf('Count: %d\n', count)
end

