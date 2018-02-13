%% load 
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
    
%%
figure; % Make a figures showing the data and the clusters
featIdx = [6, 9];

for i = unique(idx)'
    ii = find(idx == i); h = plot(+data(ii, featIdx(1)), +data(ii, featIdx(2)), 'o'); hold on;
    col = rand(1, 3); set(h, 'Color', col, 'MarkerFaceColor', col);
    xi1 = +data(i, featIdx(1)) * ones(size(ii)); xi2 = +data(i, featIdx(2)) * ones(size(ii));
    line([+data(ii, featIdx(1)), xi1]', [+data(ii, featIdx(2)), xi2]', 'Color', col);
end
axis equal tight;


%% Test lasso regression for feature selection

close all; clear all;
load 'D:\Data\Development\Projects\PhD GeoInformatics\Data\Feature Selection\Hyperspectral\BotswanaPr.mat';
data = remclass(data);
data = setprior(data, 0);
cs = classsizes(data);
myPreferredFeatures

%%
load('D:\Data\Development\Projects\PhD GeoInformatics\Data\Feature Selection\DataAllWin5NoBorder2.mat')
data = changelablist(dataAll, 'Default');
data = setprior(data, 0);
fl = cellstr(getfeatlab(data));
idx = strmatch('Lbp', fl);
data(:, idx)=[];
fl = strrep(fl, 'Ndvi', 'NDVI');
fl = strrep(fl, 'irRat', 'RVI');
fl = strrep(fl, 'IrRat', 'RVI');

%%
nl = getnlab(data);
d = prdataset(+data);
d = d*scalem(d, 'variance');
d = setlabtype(d, 'targets');
d = settargets(d, nl)

if false
    w = lassor(d, 100);

    % NB the first coeff in w is the constant/offset and must be discarded
    +w
    sum(+w > 1e-6)

    fl = cellstr(getfeatlab(dataAll));
    fl(find(+w>0)+1)
end

w = FeatSelLassoM(data, 100, 0)

fl(+w)

%% test MVFS

close all hidden; clear all;
load('D:\Data\Development\Projects\PhD GeoInformatics\Data\Feature Selection\DataAllWin5NoBorder2.mat')
dataAll = changelablist(dataAll, 'Default');
dataAll = setprior(dataAll, 0);
fl = cellstr(getfeatlab(dataAll));
idx = strmatch('Lbp', fl);
dataAll(:, idx)=[];
fl = cellstr(getfeatlab(dataAll));

cs = classsizes(dataAll);
if true
    cs(1) = cs(2);
else
    cs = min(cs)*ones(1,3);
end
randreset;
data = gendat(dataAll, cs);
data = changelablist(data, 'Default');
data = setprior(data, 0);
%%
load('D:\Data\Development\Projects\PhD GeoInformatics\Data\Feature Selection\Synthetic.mat');


%%
w = FeatSelMultiViewM(data, 5, 7)

%% Making h for MVFS
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
