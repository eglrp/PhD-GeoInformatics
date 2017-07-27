%% A reworking of my original feature clustering and ranking
%--------------------------------------------------------------------------
%% View features
%--------------------------------------------------------------------------
clear all; close all;
%load('F:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\DataAllWin5NoBorder2.mat');
load('C:\Data\Development\Projects\MSc GeoInformatics\Data\Feature Selection\DataAllWin5NoBorder2.mat')
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

%% play with bootstrap fs stability eval
close all hidden
cs = classsizes(dataAll);
%cs(1) = cs(2);
cs(1:2) = cs(3);

subData = gendat(gendat(dataAll), cs); %if N spec'd no sample with repl ??
% subData = gendat(dataAll, [5000 5000 5000]); %if N spec'd no sample with repl ??
subData = setprior(subData, 0);
% w = subData*featself([], naivebc, 6);

preferredFeatures = [9 5:8 10 1:4 19:22 15:18];

% clear res*
resi = BootstrapFsEval(subData, featseli([], naivebc, 6), 'numBootStraps', 10, 'numFeatures', 6);
resf = BootstrapFsEval(subData, featself([], naivebc, 46), 'numBootStraps', 10, 'numFeatures', 6);
resfcr = BootstrapFsEval(subData, FeatSelClusterRankM([], naivebc, 6, [], 'preferredFeatures', preferredFeatures, 'clusterThresh', 0.175), 'numBootStraps', 10, 'numFeatures', 6);
% TO DO: use a better criterion than naivebc for be i.e. something that considers features in combination (and others?)
resb = BootstrapFsEval(subData, featselb([], naivebc, 0), 'numBootStraps', 10, 'numFeatures', 6);

resr = BootstrapFsEval(subData, FeatSelReliefFM([], 20, 6), 'numBootStraps', 10, 'numFeatures', 6);

resj = BootstrapFsEval(subData, FeatSelFeastM([], 'jmi', 0), 'numBootStraps', 10, 'numFeatures', 6);
resj2 = BootstrapFsEval(subData, FeatSelFeastM([], 'jmi', 6), 'numBootStraps', 10, 'numFeatures', 6);
preferredFeatures = ones(1, size(subData, 2));
preferredFeatures([1:9]) = 1;
preferredFeatures([10 15:22]) = 2;
preferredFeatures([11:14]) = 10;
preferredFeatures([23:46]) = 10;
resj3 = BootstrapFsEval(subData(:, [1:10, 15:end]), FeatSelFeastM([], 'jmi', 6, {'preferredFeatures', ...
    preferredFeatures([1:10, 15:end])', 'similarityThresh', 0.025}), 'numBootStraps', 10, 'numFeatures', 6);

%load 'C:\Data\Development\Projects\MSc GeoInformatics\Data\Feature Selection\FeatSel2.mat'
save 'C:\Data\Development\Projects\MSc GeoInformatics\Data\FeatSel2.mat' resi resf resfcr resb resr resj

% Df = corr(resf.FeatRank(1:6,:), 'type', 'Spearman');
% pDf = triu(Df, 1);
% sum(pDf(:))/sum(1:10-1)
% Dfcr = corr(resfcr.FeatRank(1:6,:), 'type', 'Spearman');
% pDfcr = triu(Dfcr, 1);
% sum(pDfcr(:))/sum(1:10-1)
% Dr = corr(resr.FeatRank(1:6,:), 'type', 'Spearman');
% pDr = triu(Dr, 1);
% sum(pDr(:))/sum(1:10-1)
% 

% NOTES
%-----------
% - FFS & BE are ~accurate but less stable than FCR
% - Ranking is ~stable but less accurate than FCR
% - RELIEFF is similar to FCR - slightly more accurate and slightly less
% stable i.e. just use Relieff is a valid suggestion.  It will do less well
% with higher dimensional spaces though (according to the paper) and 
% doesn't allow hand picking of features.  We need to better understand how
% it does this.  

%% Investigate stability of heirarchical clustering under bootstraps
% to be able to compare eg
close all hidden
cs = classsizes(dataAll);
cs(1) = cs(2);
numBootStraps = 5;
for i = 1:numBootStraps
    fprintf('Boot strap %i of %i\n', i, numBootStraps);
    subData = gendat(gendat(dataAll), cs); %if N spec'd no sample with repl ??
    subData = setprior(subData, 0);
    res(i) = FeatureClusterRank(subData, 'clusterThresh', 0.175);
end

% compare clusters across bootstraps

fprintf('Clusters ranked by accuracy for each bootstrap\n');
for i = 1:numBootStraps
    [clustAcc_ clustIdx] = sort(res(i).ClustAcc);
    fprintf('%d ', clustIdx');
    fprintf('\n');
end

fprintf('\nCluster components for each bootstrap\n');
for ci = 1:length(res(1).ClustAcc)
    fprintf('Cluster %d\n', ci);
    for i = 1:numBootStraps
        fprintf('%d ', res(i).ClustFeatNLab{ci});
        fprintf('\n');
    end
end
% NOTE
% - The clusters sorted by number are v stable across bootstraps but not 100%
% - The clusters sorted by accuracy are mildly stable across bootstraps

%% Experiment with ReliefF algorithm and FEAST toolbox
subData2 = subData(1:5:end, :);
if true
    %ranked = feast('mrmr',6,+subData,getnlab(subData));
%     ranked = feast('jmi',6,+subData,getnlab(subData));
    % ranked = feast('fcbf',10,+subData,getnlab(subData),0.005);
    [ranked, weight] = relieff(+subData, cellstr(getlab(subData)), 20, 'prior', 'uniform', 'method', 'classification', 'categoricalx', 'off');
else
    w = FeatSelReliefFM([], 20, 6);
    [w, r] = subData2*w;
    [ranked, weight] = relieff(+subData2, cellstr(getlab(subData2)), 20, 'prior', 'uniform', 'method', 'classification', 'categoricalx', 'off');
end
% "ranked" is an index not a weight

% [weight_, sortIdx] = sort(-weight);
% fl = cellstr(getfeatlab(subData2));
fl(ranked)

clfrs = {opencvsvc([], [], {'SVMType', 'C_SVC', 'KernelType', 'RBF', 'Gamma', 2, 'C', 1, 'ClassWeights', ones(1, getsize(subData2, 3))}), ...
        opencvrtreec([], [], {'Priors', ones(1, getsize(subData2, 3))/getsize(subData2, 3), 'MaxNumOfTreesInTheForest', 5, 'NActiveVars', 4, 'CalcVarImportance', true, 'MaxDepth', 10, 'ForestAccuracy', 0.025}), ...
        opencvdtreec([], 12, {'Priors', ones(1, getsize(subData2, 3))/getsize(subData2, 3), 'MaxDepth', 12, 'Use1seRule', false, 'UseSurrogates', false, 'CVFolds', 5, 'TruncatePrunedTree', true, 'MinSampleCount', min(classsizes(subData2))/100}), ...
        naivebc};

for ci = 1:length(clfrs)
    fprintf('    Clfr eval %i of %i\n', ci, length(clfrs));

    [err, cerr, nlabOut] = prcrossval(subData(:, ranked(1:6)), scalem([], 'variance')*clfrs{ci}, 5, 1);
    c = confmat(getnlab(subData), nlabOut);
    clfCn{ci} = c./repmat(sum(c, 2), 1, size(c, 2));
    clfAcc(ci) = mean(diag(clfCn{ci}));
    clfName{ci} = struct(clfrs{ci}).name;
end
% NOTE - this produced good clf accuracies in some cases superior to FCR&FS

%% Make prtools dataset from NIPS data
read_data
%%
clear all; close all;
baseDir = 'C:\Data\Development\Projects\MSc GeoInformatics\Data\Feature Selection\NIPS\Data\';
dataDirs = {'madelon', 'gisette'};
%%
for i = 1:length(dataDirs)
    load([baseDir dataDirs{i}]);
    data = prdataset(full([X_train;X_valid]), full([Y_train;Y_valid]));
    data = setprior(data, 0);
    data = setfeatlab(data, cellstr(num2str([1:size(data, 2)]')));
    save([baseDir dataDirs{i} '_dataset.mat'], 'data');
end

%% Test FS on NIPS data
for i = 1:length(dataDirs)
    load([baseDir dataDirs{i} '_dataset.mat']);
    resj{i} = BootstrapFsEval(data, FeatSelFeastM([], 'jmi', 20), 'numBootStraps', 5, 'numFeatures', 20);
    % resr{i} = BootstrapFsEval(data, FeatSelReliefFM([], 10, 20), 'numBootStraps', 5, 'numFeatures', 20);
end    
%%
ranked = feast('jmi', 20, +data(1:5000,:), getnlab(data(1:5000,:)));

w = FeatSelClusterRankM([], naivebc, 20, [], 'clusterThresh', 0.2);
w = data*w
% resfcr = BootstrapFsEval(data, FeatSelClusterRankM([], naivebc, 10, [], 'clusterThresh', 0.175, 'showFigures', true), 'numBootStraps', 5, 'numFeatures', 10);
clfrs = {opencvsvc([], [], {'SVMType', 'C_SVC', 'KernelType', 'RBF', 'Gamma', 2, 'C', 1, 'ClassWeights', ones(1, getsize(data, 3))}), ...
        opencvrtreec([], [], {'Priors', ones(1, getsize(data, 3))/getsize(data, 3), 'MaxNumOfTreesInTheForest', 5, 'NActiveVars', 4, 'CalcVarImportance', true, 'MaxDepth', 10, 'ForestAccuracy', 0.025}), ...
        opencvdtreec([], 12, {'Priors', ones(1, getsize(data, 3))/getsize(data, 3), 'MaxDepth', 12, 'Use1seRule', false, 'UseSurrogates', false, 'CVFolds', 5, 'TruncatePrunedTree', true, 'MinSampleCount', min(classsizes(data))/100}), ...
        naivebc};

for ci = 1:length(clfrs)
    fprintf('    Clfr eval %i of %i\n', ci, length(clfrs));

    [err, cerr, nlabOut] = prcrossval(data(:,+w), scalem([], 'variance')*clfrs{ci}, 5, 1);
    c = confmat(getnlab(data), nlabOut);
    clfCn{ci} = c./repmat(sum(c, 2), 1, size(c, 2));
    clfAcc(ci) = mean(diag(clfCn{ci}));
    clfName{ci} = struct(clfrs{ci}).name;
end

%%
% preferredFeatures = ones(size(subData, 2), 1);
% preferredFeatures(1:10) = 1;
% preferredFeatures(11:22) = 2;
% preferredFeatures(23:end) = 10;
fl = cellstr(getfeatlab(subData));
preferredFeatures = ones(1, size(subData, 2));
preferredFeatures([1:9]) = 1;
preferredFeatures([10 15:22]) = 2;
preferredFeatures([11:14]) = 10;
preferredFeatures([23:46]) = 10;
[ranked, fm] = JMIC(subData, 6, 'preferredFeatures', preferredFeatures(:), 'similarityThresh', 0.1);
fl(ranked)

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