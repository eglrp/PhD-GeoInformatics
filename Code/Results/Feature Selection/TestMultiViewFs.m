%% Test AP clustering
%% load data

close all hidden; clear all;
load('D:\Data\Development\Projects\PhD GeoInformatics\Data\Feature Selection\DataAllWin5NoBorder2.mat')
dataAll = changelablist(dataAll, 'Default');
dataAll = setprior(dataAll, 0);
fl = cellstr(getfeatlab(dataAll));
idx = strmatch('Lbp', fl);
dataAll(:, idx)=[];
fl = strrep(fl, 'Ndvi', 'NDVI');
fl = strrep(fl, 'irRat', 'RVI');
fl = strrep(fl, 'IrRat', 'RVI');

data = dataAll;

%% load data option 2
close all hidden; clear all;

load('D:\Data\Development\Projects\PhD GeoInformatics\Data\Feature Selection\Hyperspectral\BotswanaPr.mat')
hsNames = {'Botswana', 'KSC'};
data = remclass(data);
data = setprior(data, 0);
fl = cellstr(getfeatlab(data));

%% prep data
data = data * scalem(data, 'variance'); % scale to unit variance
S = -distm((+data)'); % -ve euclidean distance betw feats
n = size(S, 1); % num feats
tmp = triu(S, 1) + tril(S, -1); % kind of unnecessary as Sii = 0 already
pref = sum(tmp(:)) / (n * (n - 1)); % from paper
% S = tmp + diag(pref*ones(n, 1));


%% do the clustering

[idx, netsim, i, unconverged, dpsim, expref] = apcluster(S, pref);

fprintf('Number of clusters: %d\n', length(unique(idx)));
fprintf('Fitness (net similarity): %f\n', netsim);

exemplars = unique(idx);
if false
    tmpi = idx;
    tmpi(exemplars) = 1:length(exemplars);
    idx = tmpi(idx);
    exemplars = 1:length(exemplars);
end

count = 1;
for i = exemplars
    ii = find(idx == i);
    fprintf('Cluster %d, Exemplar %sf\n', count, fl{i});
    fprintf('%s, ',fl{ii});
    fprintf('\n\n');
    count = count + 1;
end    
    % for i = 1:nclust
% %     fprintf('Cluster %d, Accuracy %.3f\n', clustIdx(i), clustAcc(clustIdx(i)));
%     fprintf('%s, ',clustFl{clustIdx(i)}{:});
%     fprintf('\n');
% end
    
%% Test lasso regression for feature selection

close all; clear all;
load 'D:\Data\Development\Projects\PhD GeoInformatics\Data\Feature Selection\Hyperspectral\BotswanaPr.mat';
data = remclass(data);
data = setprior(data, 0);
cs = classsizes(data);
myPreferredFeatures

%% data option 2
load('D:\Data\Development\Projects\PhD GeoInformatics\Data\Feature Selection\Synthetic.mat')
data = setprior(data, 0);

% There are 3 groups of 5 features, (1:5) are the correct ones, then (6:10)
% =(1:5) +noise & 11:15 = 1:5 + more noise, 16 and 17 are spurious pure noise
%i.e. we should select 5 features from the 3 groups.  There should be >5
%clusters (1,6,11), (2,7,12), (3,8,13), (4,9,14), (5, 10, 15), (16), (17).

%%
nl = getnlab(data);
d = prdataset(+data);
d = 10*d*scalem(d, 'variance');
d = setlabtype(d, 'targets');
d = settargets(d, nl);

w = FeatSelLassoM(data, 5, 0)

fl(+w)

w = featself(data, 'mi', 5)

fl(+w)

% Note: forward selection does better than LASSO

%% Test MVFS with synthetic data
clear all
load('D:\Data\Development\Projects\PhD GeoInformatics\Data\Feature Selection\Synthetic.mat');
% There are 3 groups of 5 features, (1:5) are the correct ones, then (6:10)
% =(1:5) +noise & 11:15 = 1:5 + more noise, 16 and 17 are spurious pure noise
%i.e. we should select 5 features from the 3 groups.  There should be >5
%clusters (1,6,11), (2,7,12), (3,8,13), (4,9,14), (5, 10, 15), (16), (17).
fl = cellstr(getfeatlab(data));
data = gendat(data,1000)
%%
w = FeatSelMultiViewM(data, .1, 7)
fl(+w)
% here there is not enough punishment of within view coeff's so we get
% redundant feats

w = FeatSelMultiViewM(data, 10, 7)
fl(+w)
% here we get better feats with less redundancy although still pretty
% rubbish

w = FeatSelMultiViewM(data, 100, 7)
fl(+w)
% here we get a whole lot of redundancy again, so something is wrong...



%% Experiment making h for MVFS and quadprog
% 2 clusters 1:5 and 6:10
x = 1:10;  
h = [ones(5,5), zeros(5,5);ones(5,5), zeros(5,5)]
h = [ones(5,5), zeros(5,5);zeros(5,5),ones(5,5)]
x*h*x'
sum(x(1:5)).^2 + sum(x(6:end)).^2

clustIdx{1} = [1,3,2,7,9];
clustIdx{2} = [5,4,6,8,10];
h2 = zeros(10, 10);
for i = 1:length(clustIdx)
    for j = 1:length(clustIdx{i})
        h2(clustIdx{i}(j), clustIdx{i}) = 1;
    end
end

x*h2*x'
sum(x(clustIdx{1})).^2 + sum(x(clustIdx{2})).^2
