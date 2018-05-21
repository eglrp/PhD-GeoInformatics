% test apcluster and compare to hierarchical

%% manually load data from CompareFsMethods
ii = 1;
data = cdata{ii};
clusterThresh_ = clusterThresh(ii);
numFeatures(ii)

% ii=5, for pref > median(S), num clusters = num features
% ii=4, for pref = 1.35*median(S), num ap clusters << num hier. clusters (49)
% ii=6, for pref = 1.35*median(S), num ap clusters is quite sensitive to
%   pref variation, pref = median(S) is better
% ii=1, for pref = 1.35*median(S), num ap clusters is quite sensitive to
%   pref variation, actually it is quite sensitive period.

%%
fl = cellstr(getfeatlab(data));

dataNorm = (data*scalem(data, 'variance'));
S = -((+dataNorm)' * proxm2((+dataNorm)', 'correlation'));
% c = corr(+data);
% S = abs(c);
%     S = -distm((+dataNorm)'); % -ve euclidean distance betw feats
n = size(S, 1); % num feats
tmp = triu(S, 1) + tril(S, -1); % kind of unnecessary as Sii = 0 already
pref = 2.5*sum(tmp(:)) / (n * (n - 1)); % from paper but inc slightly otherwise we dont get enough features
%pref = 1.2*median(S(:));
% pref = 0.99;
pref = 1.35*median(S(:));

S = S + 1e-9 * randn(size(S, 1), size(S, 2));
[idx, netsim, i, unconverged, dpsim, expref] = apcluster(S, pref, 'maxits', 10000, 'dampfact', 0.9);
%     [idx, netsim, i, unconverged, dpsim, expref] = apcluster(S, pref);
%     if unconverged
%         print('unconverged\n')
%         S = S + 1e-9 * randn(size(S, 1), size(S, 2));
%         [idx, netsim, i, unconverged, dpsim, expref] = apcluster(S, pref, 'maxits', 1000, 'dampfact', 0.8);
if unconverged
    error('ERROR: unconverged');
end
%     end
nclust = length(unique(idx));
tmp(unique(idx)) = 1:nclust;
lab = tmp(idx); % labels 1 indexed
fprintf('Number of clusters: %d\n', length(unique(idx)));

for i = 1:nclust
    if true
        fprintf('\nCluster %d\n', i);
        disp(fl(lab==i)')
    end
end

if false
    fprintf('Fitness (net similarity): %f\n', netsim);

    exemplars = unique(idx);
    count = 1;
    for i = exemplars'
        ii = find(idx == i);
        fprintf('Cluster %d, Exemplar %sf\n', count, fl{i});
        fprintf('%s, ',fl{ii});
        fprintf('\n');
        count = count + 1;
    end            
    fprintf('\n');
end
%% compare to hierarchical clustering

c = corr(+data);
dendg = hclust(1-abs(c), 'average'); % dendrogra

%threshold correlation at 0.2
nclust = sum(dendg(2, :) > clusterThresh_); %13
lab = hclust(1-abs(c), 'average', nclust); % labels

if false
    figure;
    plotdg(dendg)
    xidx = str2double(cellstr(get(gca, 'XTickLabel')));
    set(gca, 'XTickLabel', fl(xidx));
    % rotatetl(gca, 90, 'bottom');
    hold on
    hp = plot(0:length(xidx)+1, clusterThresh_*ones(1,length(xidx)+2), 'k--');
    axis tight
    set(gca, 'TickDir', 'out')
    set(gca, 'TickLength', [0 0])
    set(gca, 'box', 'on')
    set(gca, 'YGrid', 'on')
    set(gca, 'FontSize', 9)
    p = get(gca, 'Position');
    p(2) = p(2)+p(end)*.15;
    p(end) = p(end)*.85;
    % set(gca, 'Position', p);
    yl = get(gca, 'YTickLabel');
    yln = str2double(cellstr(yl));
    % set(gca, 'YTickLabel', max(yln)-yln);
    ylabel('Dissimilarity', 'FontSize', 11)
%         xlabel('Features', 'FontSize', 11)
%     view(90, 90)
    hl = legend(hp, 'Dissimilarity threshold', 'Location', 'NorthWest')
    fontsize(11)
    set(hl,'FontSize',11);
    view(90,90)   % for paper
end
for i = 1:nclust
    if true
        fprintf('\nCluster %d\n', i);
        disp(fl(lab==i)')
    end
end

