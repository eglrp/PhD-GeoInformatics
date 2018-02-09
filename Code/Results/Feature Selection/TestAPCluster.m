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

