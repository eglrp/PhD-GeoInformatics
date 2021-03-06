function res = FeatureClusterRank(data, varargin)

clusterMethod = 'ap';  %'ap'
apclusterCrit = 'correlation';  % proxm distance criterion for AP clustering 
apclusterParam = [];  % proxm distance criterion for AP clustering 
useExemplars = true;

clusterThresh = 0.2; %for hierarchical only
criterion = naivebc;
preferredFeatures = []; % preferred features to select from clusters in order of preference
showFigures = false;
useCorrelation = true;
clusterType = 'average';
jmiFormulation = false;  % use MI to rank features in clusters and JMI to choose clusters
fsFormulation = false;  % use FS to rank clusters

ModifyDefaultArgs(varargin);

%%
% if jmiFormulation && ~strcmpi(criterion, 'mi')

if isempty(getfeatlab(data))
    fl = num2cell(1:size(data, 2));
else  % TO DO - save this in the data set or do it outside of this call
    fl = cellstr(getfeatlab(data));
    fl = strrep(fl, 'Ndvi', 'NDVI');
    fl = strrep(fl, 'irRat', 'RVI');
    fl = strrep(fl, 'IrRat', 'RVI');
end

if strcmpi(clusterMethod, 'ap')  % affinity propagation
    fprintf('Affinity propagation:\n')
    % Note that the scaling below makes a huge difference to clustering
    % performance.  Scaling to unit variance works much better, but the
    % paper only talks of scaling to range [0,1].  BUT this scaling is not
    % good for spectral data where relative scale matters.
    % Their dist similarity measure
    % is a weakness as closely spaced features are not the only redundant
    % ones, features that are linearly scaled versions of each will be
    % redundant too but they are not close.  Then there are non-lin dependencies too.

    dataNorm = (data*scalem(data, 'variance'));
    %dataNorm = 10*(data*scalem(data, 'domain'));
    if true
        S = -((+dataNorm)' * proxm2((+dataNorm)', apclusterCrit, apclusterParam));
    %     S = -distm((+dataNorm)'); % -ve euclidean distance betw feats
        n = size(S, 1); % num feats
        tmp = triu(S, 1) + tril(S, -1); % kind of unnecessary as Sii = 0 already
        pref = 2.5*sum(tmp(:)) / (n * (n - 1)); % from paper but inc slightly otherwise we dont get enough features
        pref = median(S(:));
        [idx, netsim, i, unconverged, dpsim, expref] = apcluster(S, pref, 'maxits', 100000, 'dampfact', 0.7);
        if unconverged
            fprintf('unconverged\n')
%             S = S + 1e-9 * randn(size(S, 1), size(S, 2));
%             pref = pref*1.2;
            [idx, netsim, i, unconverged, dpsim, expref] = apcluster(S, pref, 'maxits', 100000, 'dampfact', 0.5);
            if unconverged
                error('ERROR: unconverged');
            end
        end
    elseif false
%         c = 1-abs(corr(+dataNorm));
        S = ((+dataNorm)' * proxm2((+dataNorm)', apclusterCrit, apclusterParam));
        lab = dcluste(S);
        idx = lab(:, ceil(size(lab,2)/2+1));
    else
        S = ((+dataNorm)' * proxm2((+dataNorm)', apclusterCrit, apclusterParam));
%         c = 1-abs(corr(+dataNorm));
        idx = exemplar(S, 0.2, 0.5);        
    end

    exemplars = unique(idx);
    nclust = length(exemplars);
    tmp(exemplars) = 1:nclust;
    lab = tmp(idx); % labels 1 indexed
    fprintf('Number of clusters: %d\n', length(exemplars));
    
    if showFigures
%         fprintf('Fitness (net similarity): %f\n', netsim);

%         exemplars = unique(idx);
        count = 1;
        for i = exemplars'
            ii = find(idx == i);
            fprintf('Cluster %d, Exemplar %s\n', count, fl{i});
            fprintf('%s, ',fl{ii});
            fprintf('\n');
            count = count + 1;
        end            
        fprintf('\n');
    end
else
    if useCorrelation
        c = corr(+data);
    else
        dataNorm = (10*(data*scalem(data, 'domain')));
        c = zeros(size(dataNorm, 2), size(dataNorm, 2));
%         for i = 1:size(dataNorm, 2)
%             h(i) = entropy(+dataNorm(:, i));
%         end
        for i = 1:size(dataNorm, 2)
            for j = i:size(dataNorm, 2)
                c(i, j) = kernelmi(+dataNorm(:, i), +dataNorm(:, j));
%                 c(i, j) = c(i, j)/min(h(i), h(j)); % this should have a max of 1
            end
        end
        c = c + triu(c, 1)';   % + diag(ones(1, size(data, 2)));
        c(h == 0, :) = 0;
        c(:, h == 0) = 0;
        
%     else
%         dataNorm = (10*(data*scalem(data, 'domain')));
%         c = zeros(size(dataNorm, 2), size(dataNorm, 2));
%         for i = 1:size(dataNorm, 2)
%             h(i) = entropy(+dataNorm(:, i));
%         end
%         for i = 1:size(dataNorm, 2)
%             for j = i:size(dataNorm, 2)
%                 c(i, j) = mi(+dataNorm(:, i), +dataNorm(:, j));
%                 c(i, j) = c(i, j)/min(h(i), h(j)); % this should have a max of 1
%             end
%         end
%         c = c + triu(c, 1)';   % + diag(ones(1, size(data, 2)));
%         c(h == 0, :) = 0;
%         c(:, h == 0) = 0;
    end

    if showFigures
        figure;
        hp = (imagesc(sqrt(abs(c))));
        set(gca,'YTick',1:size(c,1))
        set(gca,'YTickLabel',fl);
        set(gca,'XTick',1:size(c,2))
        set(gca,'XTickLabel',fl);
    %     rotatetl(gca,90,'top');
        colormap gray
        axis square
    end
    %heirarchical clustering
    fprintf('Heirarchical Clustering\n')

    if useCorrelation
        dendg = hclust(1-abs(c), clusterType); % dendrogra
    else
        dendg = hclust(1-abs(c), clusterType); % dendrogra
    end
    %threshold correlation at 0.2
    nclust = sum(dendg(2, :) > clusterThresh); %13
    lab = hclust(1-abs(c), clusterType, nclust); % labels
    fprintf('Number of clusters: %d\n', nclust);

    if showFigures
        figure;
        plotdg(dendg)
        xidx = str2double(cellstr(get(gca, 'XTickLabel')));
        set(gca, 'XTickLabel', fl(xidx));
        % rotatetl(gca, 90, 'bottom');
        hold on
        hp = plot(0:length(xidx)+1, clusterThresh*ones(1,length(xidx)+2), 'k--');
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
        %title(hl, '')
        legend(gca, 'boxoff')
        fontsize(11)
        set(hl,'FontSize',11);
        view(90,90)   % for paper
%         grid off
%         print("C:\Data\Development\Projects\PhD GeoInformatics\Docs\My Docs\Thesis\VHR Spekboom Canopy Cover Mapping\Fig. 6  Clustering of correlated features.png",...
%             '-dpng', '-r600')
        
    end
end % hierarchical clustering

% fprintf('K-Means Clustering')
% lab = kmeans((dataAll')*scalem(dataAll', 'variance'), nclust);
% for i = 1:nclust
%     fprintf('\nCluster %d\n', i);
%     disp(fl(lab==i))
% end

%% Get a ranking of inividual features
fprintf('\n');
if ~ismapping(criterion)
    if strcmpi(criterion, 'mi') || strcmpi(criterion, 'nmi') || strcmpi(criterion, 'jmi')
        dataNorm = (10*data*scalem(data, 'domain'));
    else
        dataNorm = data*scalem(data, 'variance');
    end
end
% clear featAcc featAccComb w
fprintf('Ranking individual features:\n-----------------------------\n');
if useExemplars && exist('exemplars', 'var')
    nf = length(exemplars);
    fi = exemplars;
else
    nf = length(fl);
    fi = 1:nf;
end
featAcc = inf*ones(length(fl), 1);
for i = 1:nf
    fprintf('%d,',fi(i));
    try
        if ~ismapping(criterion) && strcmpi(criterion, 'distcorr')
            featAcc(fi(i)) = -FaDCor(+dataNorm(:, fi(i)), getnlab(dataNorm));
        elseif ~ismapping(criterion)
            featAcc(fi(i)) = -FeatEvalMi(dataNorm(:, fi(i)), criterion);
        else
            %featAcc is error rate
            [featAcc(fi(i)), cerr, nlabOut, stds, r] = prcrossval(data(:, fi(i)), criterion, 5);
        end
    catch
        fprintf('Error');
        featAcc(fi(i)) = inf;
    end
end
fprintf('\n');

[featAcc_, featIdx] = sort(featAcc);
% fl(featIdx)
res.FeatIndividualAcc = featAcc;

%% Combine feature ranking with clustering

fprintf('Feature ranking and clustering:\n-----------------------------\n');
if showFigures
    fprintf('Name\tCluster\tAcc\n');
    clear table;
    table(1, :) = {'Name','Cluster', 'Acc'};
    for i = 1:length(featAcc)
        table(i+1,:) = {fl{featIdx(i)}, lab(featIdx(i)), featAcc(featIdx(i))};
        %fprintf('%s\t%d\t%.2f\n', fl{featIdx(i)}, lab(featIdx(i)), featAcc(featIdx(i)));
    end
    disp(table);
end

%find average accuracy per cluster
clustAcc = [];
for i = 1:nclust
    if useExemplars && exist('exemplars', 'var')
        clustAcc(i) = featAcc(exemplars(i));  %min(featAcc(lab==i));  %
    else
        clustAcc(i) = median(featAcc(lab==i));
    end
end

[clustAcc_ clustIdx] = sort(clustAcc);

%renumber the clusters by score to try and make them consistent between
%bootstraps (actually I doubt this will improve things any more than
%the output clusters from hclust which are ordered by correlation somehow
%but it is nice to have the clusters in order of their scores)
if true
    clustRank(clustIdx) = 1:nclust;
    lab = clustRank(lab);  % equivalent?
    res.FeatClusterNLab = lab; % cluster numbers for each feature
    %res.ClustFeatNLab = res.ClustFeatNLab{clustIdx};
    if useExemplars && exist('exemplars', 'var')
        exemplars = exemplars(clustIdx);
    end
    clustIdx_ = clustIdx;
    clustIdx = 1:nclust;  % 
    res.ClustAcc = clustAcc_;
else
    res.ClustAcc = clustAcc;
end

clustFl = {};
clustFeatNLab = {};
for i = 1:nclust
    if showFigures
        fprintf('\nCluster %d\n', i);
        disp(fl(lab==i))
    end
    clustFl(i) = {fl(lab==i)};
    clustFeatNLab{i} = find(lab==i);
end
res.ClustFeatNLab = clustFeatNLab;  % feature numbers in each cluster
res.ClustFeatLab = clustFl;         % feature labels in each cluster

bestFeats = zeros(1, length(res.ClustAcc));
% find the most accurate feature from each cluster
for i = 1:length(res.ClustAcc)
    clusterFeatIdx = res.ClustFeatNLab{i}; %clustIdx(i)};
    for j = 1:length(clusterFeatIdx)  % check if there are any "preferred" features in this cluster
        fi = find(clusterFeatIdx(j) == preferredFeatures, 1);
        if ~isempty(fi)
            bestFeats(i) = preferredFeatures(fi(1)); % choose the feature listed first
            break;
        end
    end
    if bestFeats(i) == 0  % else use the "best" (highest scored) feature
        if useExemplars && exist('exemplars', 'var')
            bestFeats(i) = exemplars(i);
        else
            [featAcc_ sortClusterFeatIdx] = sort(res.FeatIndividualAcc(clusterFeatIdx));
            bestFeats(i) = clusterFeatIdx(sortClusterFeatIdx(1));
        end
    end
end

if jmiFormulation    % select clusters using JMI rather than plain ranking
    dataNorm = (10*(data*scalem(data, 'domain')));
    jmiFeats = feast('jmi', length(bestFeats), +dataNorm(:, bestFeats), getnlab(dataNorm));
    clustIdx = jmiFeats;
end

if fsFormulation    % select clusters using forward selection with the spec'd criterion
    if false
        dataNorm = (10*(data*scalem(data, 'variance')));
        w = featself(dataNorm(:, bestFeats), opencvknnc([], 3), length(bestFeats));
    else
        if ismapping(criterion)
            dataNorm = (10*(data*scalem(data, 'variance')));
            w = featself(dataNorm(:, bestFeats), criterion, length(bestFeats), 5);
        else
            dataNorm = (10*(data*scalem(data, 'domain')));
            w = featself(dataNorm(:, bestFeats), criterion, length(bestFeats));
        end
    end
    clustIdx = +w;
end

if showFigures
    fprintf('Cluster Ranking:\n----------------\n');
    for i = 1:nclust
        fprintf('Cluster %d, Accuracy %.3f\n', clustIdx(i), 100 - 100*clustAcc_(clustIdx(i)));
        fprintf('\t%s, ',clustFl{clustIdx(i)}{:});
        fprintf('\n');
    end
end
% ranks and scores for measuring stability
res.FeatClustScore = res.ClustAcc(res.FeatClusterNLab);


if false
    clustRank(clustIdx) = 1:length(clustIdx);
else % fractional ranking
    ttlFeats = 0;
    for i = 1:length(clustIdx)
        clustFeats = clustFeatNLab{clustIdx(i)};
        clustRank(clustIdx(i)) = (length(clustFeats)+1)/2 + ttlFeats;
        ttlFeats = ttlFeats + length(clustFeats);
    end
end

res.FeatClustRank = clustRank(res.FeatClusterNLab);
res.FeatClustIdx = clustIdx(res.FeatClusterNLab);

% [clustAcc_ clustIdx] = sort(res.ClustAcc);
% 
% fprintf('Cluster Ranking:\n----------------\n');
% 
% feats = zeros(1, length(res.ClustAcc));
% % choose the most accurate feature from each cluster
% for i = 1:length(res.ClustAcc)
%     clusterFeatIdx = res.ClustFeatNLab{clustIdx(i)};
%     for j = 1:length(clusterFeatIdx)  % check if there are any "preferred" features in this cluster
%         fi = find(clusterFeatIdx(j) == preferredFeatures, 1);
%         if ~isempty(fi)
%             feats(i) = preferredFeatures(fi(1));
%             break;
%         end
%     end
%     if feats(i) == 0  % else use the "best" (highest scored) feature
%         [featAcc_ sortClusterFeatIdx] = sort(res.FeatIndividualAcc(clusterFeatIdx));
%         feats(i) = clusterFeatIdx(sortClusterFeatIdx(1));
%     end
% end

res.Feats = bestFeats(clustIdx);
% for i = 1:nclust
% %     fprintf('Cluster %d, Accuracy %.3f\n', clustIdx(i), clustAcc(clustIdx(i)));
%     fprintf('%s, ',clustFl{clustIdx(i)}{:});
%     fprintf('\n');
% end