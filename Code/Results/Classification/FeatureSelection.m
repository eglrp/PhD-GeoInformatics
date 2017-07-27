%% Feature selection
%--------------------------------------------------------------------------
%% View features
%--------------------------------------------------------------------------
clear all; close all;
load('F:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\DataAllWin5NoBorder2.mat');
dataAll = changelablist(dataAll, 'Default');
fl = cellstr(getfeatlab(dataAll));
%remove LBP
idx = strmatch('Lbp', fl);
dataAll(:,idx)=[];

figure;
scatterd(dataAll(:, [2 5]), 'legend')

%look particularly at separation betw sb (res) and bg (blue)
scatterdui(dataAll)

%NOTES on visual inspection
%--------------------------------------------------------------------------
%- rN and NDVI correlated
%- bN is good at separating tree from sb
%- gN separated bg into a number of clusters
%- NDVI and nirN are correlated
%- NDVI and irRat are non-lin correlated (exact relationship)
%- tc1-2 useful
%- tc4 perhaps good at separating sb and tree (as for bN)
%- pc distr similar to tc
%- rc1 correlated with NDVI
%- EntropyPc1 separated tree from other quite well but has weird
%quantisation issues
%- StdPc1, KurtosisPc1 & SkewnessPc1 do not look useful
%- LbpPc1 weirdlyu quantised and seemingly not.  Other lbp features hard to
%judge because of quantisation
%- NDVI and meanNdvi correlated
%- Texture feats other than Pc1 don't appear useful 
%- RGBIR are fairly normally distr, rgG are weirdly distr

%CONCL on visual inspection
%-------------------------------------------------------------------------
%- Investigate quantisation
%- rc*, irRat, Std*, Kurtosis*, Skewness*, Lbp??? can be excl.  All texture
%feats other than *Pc1 can be excl


%% Do a correlation/cluster analysis of feature groups (as opposed to objects which would be the norm)
%%NNB change line below to dataAll
close all hidden
cs = classsizes(dataAll);
cs(1) = cs(2);

subData = gendat(gendat(dataAll), cs); %if N spec'd no sample with repl ??
% subData = gendat(dataAll, [5000 5000 5000]); %if N spec'd no sample with repl ??
subData = setprior(subData, 0);


fl = cellstr(getfeatlab(dataAll));

subData = gendat(dataAll, cs);
subData = setprior(subData, 0);
% FeatureClusterRankBe(subData);
FeatureClusterRank(subData, 'clusterThresh', 0.175);

% Cluster 7, Accuracy 0.003
% 	rN, 	nirN, 	NDVI, 	RVI, 	tc2, 	pc2, 	rc1, 	MeanRVI, 	MedianRVI, 	MeanNDVI, 	MedianNDVI, 
% Cluster 10, Accuracy 0.002
% 	StdPc1, 
% Cluster 6, Accuracy 0.003
% 	pc4, 
% Cluster 4, Accuracy 0.004
% 	rc2, 	rc4, 
% Cluster 1, Accuracy 0.017
% 	R, 	G, 	B, 	NIR, 	tc1, 	pc1, 	MeanPc1, 	MedianPc1, 
% Cluster 2, Accuracy 0.008
% 	EntropyPc1, 
% Cluster 8, Accuracy 0.002
% 	gN, 	MeanGn, 	MedianGn, 

%%
%remove LBP
% idx = strmatch('Lbp', fl);
% dataAll(:,idx)=[];

fl = cellstr(getfeatlab(dataAll));
fl = strrep(fl, 'Ndvi', 'NDVI');
fl = strrep(fl, 'irRat', 'RVI');
fl = strrep(fl, 'IrRat', 'RVI');

c = corr(+subData);

figure;
hp = imagesc(1-abs(c));
set(gca,'YTick',1:size(c,1))
set(gca,'YTickLabel',fl);
set(gca,'XTick',1:size(c,2))
set(gca,'XTickLabel',fl);
rotatetl(gca,90,'top');
colormap gray
axis square

%heirarchical clustering
dendg = hclust(1-abs(c), 'average'); % dendrogra

%threshold correlation at 0.2
thr = 0.2;
nclust = sum(dendg(2,:)>thr); %13
lab = hclust(1-abs(c), 'average', nclust); % labels

figure;
plotdg(dendg)
xidx = str2double(cellstr(get(gca, 'XTickLabel')));
set(gca, 'XTickLabel', fl(xidx));
% rotatetl(gca, 90, 'bottom');
hold on
plot(0:length(xidx)+1, thr*ones(1,length(xidx)+2), 'r-');
axis tight
set(gca, 'TickDir', 'out')
set(gca, 'TickLength', [0 0])
set(gca, 'box', 'on')
set(gca, 'YGrid', 'on')
p = get(gca, 'Position');
p(2) = p(2)+p(end)*.15;
p(end) = p(end)*.85;
set(gca, 'Position', p);
yl = get(gca, 'YTickLabel');
set(gca, 'YTickLabel', 1-str2double(cellstr(yl)));
ylabel('|Correlation Coefficient|', 'FontSize', 10)


figure;
plotdg(dendg)
xidx = str2double(cellstr(get(gca, 'XTickLabel')));
set(gca, 'XTickLabel', fl(xidx));
% rotatetl(gca, 90, 'bottom');
hold on
plot(0:length(xidx)+1, thr*ones(1,length(xidx)+2), 'r-');
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
    fprintf('\nCluster %d\n', i);
    disp(fl(lab==i))
    clustFl(i) = {fl(lab==i)};
end

% fprintf('K-Means Clustering')
% lab = kmeans((dataAll')*scalem(dataAll', 'variance'), nclust);
% for i = 1:nclust
%     fprintf('\nCluster %d\n', i);
%     disp(fl(lab==i))
% end

%% Get a ranking of inividual features
subData = setprior(subData, 0);
[tr ts] = gendat(subData, 0.5);
ktr = gendat(tr, [100 100 100]);
% fl = cellstr(getfeatlab(dataAll));

fl = cellstr(getfeatlab(subData));
fl = strrep(fl, 'Ndvi', 'NDVI');
fl = strrep(fl, 'irRat', 'RVI');
fl = strrep(fl, 'IrRat', 'RVI');

fprintf('\n');
clear featAcc featAccComb w
for i = 1:length(fl)
    fprintf('%d,',i);
    try
        if false
            featAcc(i, :) = 1-feateval(tr(:, i), 'in-in');
            featAccComb(i,1) = featAcc(i, :);
%             featAcc(i, :) = ts(:, i)*w*testc;
%             featAccComb(i,1) = (ts(:, i)*[w{:}]*maxc)*testc;
        elseif true
            %w{1} = ktr(:, i)*knnc([], 1);
             w{1} = tr(:, i)*naivebc;
%             w{1} = tr(:, i)*qdc;
            featAcc(i, :) = ts(:, i)*w*testc;
            featAccComb(i,1) = (ts(:, i)*[w{:}]*maxc)*testc;
        elseif false
            w{1} = ktr(:, i)*knnc([], 1);
            w{2} = tr(:, i)*naivebc;
    %         w{3} = tr(:, i)*fisherc;
            w{3} = tr(:, i)*qdc;
            w{4} = tr(:, i)*nmc;
            featAcc(i, :) = ts(:, i)*w*testc;
            featAccComb(i,1) = (ts(:, i)*[w{:}]*maxc)*testc;
        end
    catch
        fprintf('Error');        
        featAcc(i) = 1;
    end
end
fprintf('\n');

if 0
ktr = gendat(tr, [5000 5000 5000]);
% ktr = setprior()
wrf = librandomforestc(ktr, 50);
% wrf = librandomforestc(tr, 100, 4)
tmp = +wrf;
featAcc(:, 5) = 0.5 - tmp.importance(:,4);
end
%NB HACK use rand forest ranking
% featAccComb = featAcc(:, 5);

%sort according to ???
[featAcc_ featIdx] = sort(featAccComb);
fl(featIdx)

figure
plot(featAcc(featIdx,:))
hold all
plot(featAccComb(featIdx,:))
legend({'nn','naivebc','qdc','nmc','librandomforestc','maxc'})

set(gca,'XTick',1:length(fl))
set(gca,'XTickLabel',fl(featIdx));
rotatetl(gca,90,'bottom');

%% Combine feature ranking with clustering

fprintf('Feature ranking and clusterinF:\n-----------------------------\n');
fprintf('Name\tCluster\tAcc\n');    
clear table;
table(1,:) = {'Name','Cluster', 'Acc'};
for i = 1:length(featAcc_)
    table(i+1,:) = {fl{featIdx(i)}, lab(featIdx(i)), featAccComb(featIdx(i))};
    fprintf('%s\t%d\t%.2f\n', fl{featIdx(i)}, lab(featIdx(i)), featAccComb(featIdx(i)));
end
disp(table);

%find average accuracy per cluster
clustAcc = [];
for i = 1:nclust
    clustAcc(i) = median(featAccComb(lab==i));
end

[clustAcc_ clustIdx] = sort(clustAcc);

fprintf('Cluster RankinF:\n----------------\n');
for i = 1:nclust
    fprintf('Cluster %d, Accuracy %.3f\n', clustIdx(i), clustAcc(clustIdx(i)));
    fprintf('\t%s, ',clustFl{clustIdx(i)}{:});
    fprintf('\n');
end
for i = 1:nclust
%     fprintf('Cluster %d, Accuracy %.3f\n', clustIdx(i), clustAcc(clustIdx(i)));
    fprintf('%s, ',clustFl{clustIdx(i)}{:});
    fprintf('\n');
end

%% What sort of acc do we get if we choose one feature from each of the 1st N clusters
load('F:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\DataAllWin5NoBorder2.mat')
dataAll = changelablist(dataAll, 'Default')
dataAll = setprior(dataAll, 0);
%remove LBP
% idx = strmatch('Lbp', fl);
% dataAll(:,idx)=[];

global feats
feats = [9 15 20 23 7 6]; %ranked cluster
% wfso = featself(subData, opencvdtreec([], 12, {'Priors', [0.25 0.5 0.25], 'MaxDepth', 12, 'Use1seRule', false, ...
%         'UseSurrogates', false, 'CVFolds', 5, 'TruncatePrunedTree', true, 'MinSampleCount', min(classsizes(subData))/100}), 6);
% feats = [16    26    14     7    38    19]
wfso = featself(subData, naivebc, 6);
feats = [16    14     7    23    45    21]
wro = featrank(subData, naivebc)
feats = wro(1:6)
wbeo = featselb(subData, naivebc, 6);
feats = [2    19    21    23    31    32]
% wbeo = featselb(subData, naivebc, 6);
% feats = [16 23 21 25 38 20]; %wfso = featself(subData, libsvc([], proxm([], 'r', 2), 20), 6);
% feats = +wfso
    'NDVI'
    'pc1'
    'rc2'
    'EntropyPc1'
    'bN'
    'gN'
% feats = [9 15 23 7 6];
subData = gendat(dataAll, [2000 2000 2000]);
if false % 2 class
    nlab = getnlab(subData);
    otherIdx = getnlab(subData) == getclassi(subData, 'Background') | getnlab(subData) == getclassi(subData, 'Tree');
    ll = cellstr(getlab(subData));
    ll(otherIdx) = {'Other'};
    subData = addlabels(subData, char(ll), '2Class');
    subData = changelablist(subData, '2Class');
    subData = setprior(subData, 0);
else
    subData = changelablist(subData, 'Default');
    subData = setprior(subData, 0);
end
% [err, cerr, nlabOut] = prcrossval(subData(:, feats), opencvrtreec([], [], {'Priors', [1 3 1]/5, ...
%     'MaxNumOfTreesInTheForest', 5, 'NActiveVars', 4, 'CalcVarImportance', true, 'MaxDepth', 10, 'ForestAccuracy', 0.025}), 5, 1);
[err, cerr, nlabOut] = prcrossval(subData(:, feats), opencvdtreec([], 12, {'Priors', [0.25 0.5 0.25], 'MaxDepth', 12, 'Use1seRule', false, ...
        'UseSurrogates', false, 'CVFolds', 5, 'TruncatePrunedTree', true, 'MinSampleCount', min(classsizes(subData))/100}), 5, 1);
% [err, cerr, nlabOut] = prcrossval(subData(:, feats), opencvdtreec([], 10, {'Priors', [0.25 0.75], 'MaxDepth', 10, 'Use1seRule', false, ...
%         'UseSurrogates', false, 'CVFolds', 5, 'TruncatePrunedTree', true, 'MinSampleCount', min(classsizes(subData))/100}), 5, 1); % 2class
% [err, cerr, nlabOut] = prcrossval(subData(:, feats), scalem([], 'variance')*opencvsvc([], [], ...
%     {'SVMType', 'C_SVC', 'KernelType', 'RBF', 'Gamma', 5, 'C', 1, 'ClassWeights', double([1; 2; 1])}), 5, 1);
% [err, cerr, nlabOut] = prcrossval(subData(:, feats), scalem([], 'variance')*libsvc([], proxm([], 'r', 2), 20), 5, 1);
% [err, cerr, nlabOut] = prcrossval(subData(:, feats), scalem([], 'variance')*librandomforestc([], 100), 5, 1);
% [err, cerr, nlabOut] = prcrossval(dataAll(:, feats),qdc, 5, 1);
% [err, cerr, nlabOut] = prcrossval(subData(:, feats), dtc, 5, 1);

c = confmat(getnlab(subData), nlabOut);
cn = c./repmat(sum(c, 2), 1, size(c, 2))
1 - mean(diag(cn))

[CRn] = ReduceConfMat(cn, {[1 3], [2]}, true)
1 - mean(diag(CRn))

global w
w = scalem([], 'variance')*librandomforestc([], 100);
w = scalem([], 'variance')*libsvc([], proxm([], 'r', 2), 20);
w = opencvrtreec([], [], {'Priors', [1 3 1]/5, ...
    'MaxNumOfTreesInTheForest', 5, 'NActiveVars', 4, 'CalcVarImportance', true, 'MaxDepth', 10, 'ForestAccuracy', 0.025})
w = opencvdtreec([], 12, {'Priors', [0.25 0.5 0.25], 'MaxDepth', 12, 'Use1seRule', false, ...
        'UseSurrogates', false, 'CVFolds', 5, 'TruncatePrunedTree', true, 'MinSampleCount', min(classsizes(subData))/100})
w = opencvdtreec([], 12, {'Priors', [0.25 0.75], 'MaxDepth', 10, 'Use1seRule', false, ...
        'UseSurrogates', false, 'CVFolds', 5, 'TruncatePrunedTree', true, 'MinSampleCount', min(classsizes(subData))/100})
% w = opencvdtreec([], 12, {'Priors', [0.25 0.5 0.25], 'MaxDepth', 6, 'Use1seRule', false, ...
%         'UseSurrogates', false, 'CVFolds', 5, 'TruncatePrunedTree', true, 'MinSampleCount', min(classsizes(subData))/50})
    
w = scalem([], 'variance')*opencvsvc([], [], ...
    {'SVMType', 'C_SVC', 'KernelType', 'RBF', 'Gamma', 5, 'C', 1, 'ClassWeights', double([1; 2; 1])})
w = subData(:, feats)*w;
clfr = +w;
% clfr.save('F:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\configDTree.yaml');

clfrSvm = +w;
clfrSvm = +clfrSvm{2};
clfrSvm.save('F:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\configSvm.yaml');

% save('F:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\All_Svc_clfr3.mat', 'w', 'feats')
% save('F:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\All_RandomForest_clfr_3class.mat', 'w', 'feats')
% save('F:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\All_OpenCvRTrees_clfr_3class.mat', 'w', 'feats')
% save('F:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\All_OpenCvDTree_clfr_3class.mat', 'w', 'feats')
% save('F:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\All_OpenCvSvc_clfr_3class.mat', 'w', 'feats')
%save('F:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\All_OpenCvDTree_clfr_3class.mat', 'w', 'feats', 'subData')

s = [];
s.tr = single(+subData(:, feats));
s.trLab = single(+getnlab(subData));

%copy prev file contents before exec this line
% cv.FileStorage('F:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\trDataDTree.yaml', s);

%--------------------------------------------------------------------------
%% Eval stability of various feat selection criteria in fs
clear all; close all;
load('F:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\DataAllWin5NoBorder2.mat');
%remove LBP
fl = cellstr(getfeatlab(dataAll));
idx = strmatch('Lbp', fl);
dataAll(:,idx)=[];

fl = cellstr(getfeatlab(dataAll));
fl = strrep(fl, 'Ndvi', 'NDVI');
fl = strrep(fl, 'irRat', 'RVI');
fl = strrep(fl, 'IrRat', 'RVI');


dataAll = changelablist(dataAll, 'Default');
dataAll = setprior(dataAll, 0);
cs = classsizes(dataAll);
cs(1)=cs(2);
% subData = gendat(dataAll, [3000 3000 3000]);
subData = gendat(dataAll, cs);
subData = setprior(subData, 0);

%%
crit = {'in-in', 'maha-s', naivebc, opencvdtreec, qdc};

for i = 1:length(crit)
    i
    if (isstr(crit{i}))
        [wfs{i} rfs{i}] = featself(subData, crit{i}, 6);
    else
        [wfs{i} rfs{i}] = featself(subData, crit{i}, 6, 5);
    end
end
% fl = cellstr(getfeatlab(dataAll))
for i = 1:length(rfs)
    res(i,:) = fl(rfs{i}(:,3));
    disp(crit{i});
    fl(rfs{i}(:,3))
end

%% Eval stability of feat selection on data set sampling

for i = 1:5
    i
%     subData = gendat(dataAll, [5000 5000 5000]);
    subData = gendat(gendat(dataAll), cs);
    subData = setprior(subData, 0);
    [wfsd{i} rfsd{i}] = featself(subData, naivebc, 6, 5);
end
for i = 1:length(rfsd)
    res(i,:) = fl(rfsd{i}(:,3));
    fl(rfsd{i}(:,3))
end

%% Eval stability of feat selection on search method
% subData = gendat(dataAll, [5000 5000 5000]);
subData = setprior(subData, 0);

featSel = {@featself, @featselb, @featseli, @featselo};
for i = 1:length(featSel)-1
    [wfss{i} rfss{i}] = feval(featSel{i}, subData, naivebc, 6, 5);
end

% i = length(featSel);
% [wfss{i} rfss{i}] = feval(featSel{i}, gendat(subData, [1650 1650 1650]), opencvdtreec, 5); %full data set with cross valid too slow
d = subData;
tic
p = [1 1 1];
w = opencvdtreec([], 12, {'Priors', p./sum(p), 'MaxDepth', 12, 'Use1seRule', false, ...
        'UseSurrogates', false, 'CVFolds', 5, 'TruncatePrunedTree', true, 'CalcVarImportance', true, 'MinSampleCount', min(classsizes(subData))/100});
w = subData*w;
subData*w*testc

imp = struct(w).data{1}.getVarImportance;
[imp_ impIdx] = sort(imp, 'descend')
fl(impIdx(1:6))
% [cc ccr] = ClfrPerfMeas(d, nlabOut);

res = {};
for i = 1:length(wfss)-1
    disp(featSel{i})
    fl_ = cellstr(getfeatlab(subData*wfss{i}))
    fl_ = strrep(fl_, 'Ndvi', 'NDVI');
    fl_ = strrep(fl_, 'irRat', 'RVI');
    fl_ = strrep(fl_, 'IrRat', 'RVI');
    
    res(i,:) = fl_;
end
res(i+1,:) = fl(impIdx(1:6));
% save('D:\Data\Development\Projects\MSc GeoInformatics\Code\Results\Classification\FeatSelVariation2.mat', 'wf*', 'rf*', 'subData', 'crit', 'featSel')

%--------------------------------------------------------------------------
%%
% 
% 
% Cluster Ranking: (dist-maha)
% ----------------
% Cluster 3, Accuracy -1.199
% 	rN, 	nirN, 	NDVI, 	RVI, 	tc2, 	pc2, 	rc1, 	MeanRVI, 	MedianRVI, 	MeanNDVI, 	MedianNDVI, 
% Cluster 5, Accuracy 0.275
% 	bN, 
% Cluster 7, Accuracy 0.325
% 	rc2, 	rc4, 
% Cluster 2, Accuracy 0.412
% 	EntropyPc1, 
% Cluster 1, Accuracy 0.438
% 	R, 	G, 	B, 	NIR, 	tc1, 	pc1, 	MeanPc1, 	MedianPc1, 
% Cluster 8, Accuracy 0.505
% 	tc4, 	rc3, 
% Cluster 4, Accuracy 0.685
% 	gN, 	MeanGn, 	MedianGn, 
% Cluster 10, Accuracy 0.727
% 	EntropyRVI, 	StdRVI, 	EntropyNDVI, 	StdNDVI, 
% Cluster 9, Accuracy 0.770
% 	pc4, 
% Cluster 12, Accuracy 0.778
% 	StdPc1, 
% Cluster 11, Accuracy 0.819
% 	EntropyGn, 	StdGn, 
% Cluster 6, Accuracy 0.894
% 	tc3, 	pc3, 
% Cluster 15, Accuracy 0.903
% 	SkewnessRVI, 	SkewnessNDVI, 
% Cluster 16, Accuracy 0.991
% 	SkewnessGn, 
% Cluster 13, Accuracy 0.993
% 	SkewnessPc1, 
% Cluster 14, Accuracy 0.994
% 	KurtosisRVI, 	KurtosisNDVI, 
% Cluster 18, Accuracy 0.998
% 	KurtosisPc1, 
% Cluster 17, Accuracy 0.998
% 	KurtosisGn, 
% 
%     
% Cluster Ranking: (qdc)
% ----------------
% Cluster 3, Accuracy 0.313
% 	rN, 	nirN, 	NDVI, 	RVI, 	tc2, 	pc2, 	rc1, 	MeanRVI, 	MedianRVI, 	MeanNDVI, 	MedianNDVI, 
% Cluster 1, Accuracy 0.358
% 	R, 	G, 	B, 	NIR, 	tc1, 	pc1, 	MeanPc1, 	MedianPc1, 
% Cluster 7, Accuracy 0.391
% 	rc2, 	rc4, 
% Cluster 2, Accuracy 0.410
% 	EntropyPc1, 
% Cluster 5, Accuracy 0.424
% 	bN, 
% Cluster 4, Accuracy 0.437
% 	gN, 	MeanGn, 	MedianGn, 
% Cluster 8, Accuracy 0.479
% 	tc4, 	rc3, 
% Cluster 9, Accuracy 0.521
% 	pc4, 
% Cluster 10, Accuracy 0.526
% 	EntropyRVI, 	StdRVI, 	EntropyNDVI, 	StdNDVI, 
% Cluster 6, Accuracy 0.534
% 	tc3, 	pc3, 
% Cluster 12, Accuracy 0.538
% 	StdPc1, 
% Cluster 11, Accuracy 0.561
% 	EntropyGn, 	StdGn, 
% Cluster 15, Accuracy 0.585
% 	SkewnessRVI, 	SkewnessNDVI, 
% Cluster 16, Accuracy 0.632
% 	SkewnessGn, 
% Cluster 13, Accuracy 0.649
% 	SkewnessPc1, 
% Cluster 17, Accuracy 0.650
% 	KurtosisGn, 
% Cluster 14, Accuracy 0.654
% 	KurtosisRVI, 	KurtosisNDVI, 
% Cluster 18, Accuracy 0.656
% 	KurtosisPc1,   