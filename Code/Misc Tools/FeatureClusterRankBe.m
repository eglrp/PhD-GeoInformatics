function FeatureClusterRank(data, varargin)

clusterThresh = 0.2;
numFeats = 6;

ModifyDefaultArgs(varargin);
fl = cellstr(getfeatlab(data));
fl = strrep(fl, 'Ndvi', 'NDVI');
fl = strrep(fl, 'irRat', 'RVI');
fl = strrep(fl, 'IrRat', 'RVI');

feats = [];

c = corr(+data);

% figure;
% hp = imagesc(1-abs(c));
% set(gca,'YTick',1:size(c,1))
% set(gca,'YTickLabel',fl);
% set(gca,'XTick',1:size(c,2))
% set(gca,'XTickLabel',fl);
% rotatetl(gca,90,'top');
% colormap gray
% axis square

%heirarchical clustering
dendg = hclust(1-abs(c), 'average'); % dendrogra

%threshold correlation at 0.2
nclust = sum(dendg(2,:)>clusterThresh); %13
lab = hclust(1-abs(c), 'average', nclust); % labels

figure;
plotdg(dendg)
xidx = str2double(cellstr(get(gca, 'XTickLabel')));
set(gca, 'XTickLabel', fl(xidx));
% rotatetl(gca, 90, 'bottom');
hold on
plot(0:length(xidx)+1, clusterThresh*ones(1,length(xidx)+2), 'r-');
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
set(gca, 'YTickLabel', 1-str2double(cellstr(yl)));
ylabel('|Correlation Coefficient|', 'FontSize', 9)
view(90,90)

fprintf('Heirarchical Clustering')
clustFl = {};
for i = 1:nclust
%     fprintf('\nCluster %d\n', i);
%     disp(fl(lab==i))
    clustFl(i) = {fl(lab==i)};
end


%get pc feat for each cluster
clustFeats = prdataset(data(:,1:nclust)); %copies labels
for i = 1:nclust
%     fprintf('\nCluster %d\n', i);
%     disp(fl(lab==i))
    clustData = data(:, find(lab==i));
    w = pcam(clustData, 1);
    clustFeats(:, i) = clustData*w;
end

[featAccFull, cerr, nlabOut, stds, r] = prcrossval(clustFeats, naivebc, 5);

[tr ts] = gendat(clustFeats, 0.5);
wFull = tr*naivebc;
featAccFull = ts*wFull*testc;
%rank each cluster
for i = 1:nclust
    clustFeats_ = ts;
    clustFeats_(:,i) = clustFeats_(randperm(size(clustFeats_,1)), i);
    featAcc(i) = clustFeats_*wFull*testc;
    clustScore(i) = featAcc(i) - featAccFull;
%     [featAcc(i), cerr, nlabOut, stds, r] = prcrossval(clustFeats_, opencvdtreec, 5);
%     clustScore(i) = featAcc(i) - featAccFull;
end

[clustScore clustIdx] = sort(clustScore, 'descend');

% fprintf('K-Means Clustering')
% lab = kmeans((dataAll')*scalem(dataAll', 'variance'), nclust);
% for i = 1:nclust
%     fprintf('\nCluster %d\n', i);
%     disp(fl(lab==i))
% end


%% Combine ranking with clustering

fprintf('Feature ranking and clustering:\n-----------------------------\n');
fprintf('Name\tCluster\tAcc\n');    



fprintf('Cluster Ranking:\n----------------\n');
for i = 1:nclust
    fprintf('Cluster %d, Accuracy %.3f\n', clustIdx(i), clustScore(clustIdx(i)));
    fprintf('\t%s, ',clustFl{clustIdx(i)}{:});
    fprintf('\n');
end

% 
% for i = 1:nclust
% %     fprintf('Cluster %d, Accuracy %.3f\n', clustIdx(i), clustAcc(clustIdx(i)));
%     fprintf('%s, ',clustFl{clustIdx(i)}{:});
%     fprintf('\n');
% end