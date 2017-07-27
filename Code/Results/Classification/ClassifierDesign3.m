%% Feature selection
%--------------------------------------------------------------------------
%% View features
%--------------------------------------------------------------------------
clear all; close all;
load('G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\DataAllWin5NoBorder2.mat');
dataAll = changelablist(dataAll, 'Default');

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
fl = cellstr(getfeatlab(dataAll));
%remove LBP
idx = strmatch('Lbp', fl);
dataAll(:,idx)=[];

fl = cellstr(getfeatlab(dataAll));
c = corr(+dataAll);

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
nclust = sum(dendg(2,:)>0.2); %13
lab = hclust(1-abs(c), 'average', nclust); % labels

figure;
plotdg(dendg)
xidx = str2double(cellstr(get(gca, 'XTickLabel')));
set(gca, 'XTickLabel', fl(xidx));
rotatetl(gca, 90, 'bottom');
hold on
plot(xidx', 0.2*ones(1,length(xidx)), 'r-');
axis tight

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
dataAll = setprior(dataAll, 0);
[tr ts] = gendat(dataAll, 0.5);
ktr = gendat(tr, [100 100 100]);
fl = cellstr(getfeatlab(dataAll));

fprintf('\n');
clear featAcc featAccComb w
for i = 1:length(fl)
    fprintf('%d,',i);
    try
        w{1} = ktr(:, i)*knnc([], 1);
        w{2} = tr(:, i)*naivebc;
%         w{3} = tr(:, i)*fisherc;
        w{3} = tr(:, i)*qdc;
        w{4} = tr(:, i)*nmc;
        featAcc(i, :) = ts(:, i)*w*testc;
        featAccComb(i,1) = (ts(:, i)*[w{:}]*maxc)*testc;
    catch
        fprintf('Error');        
        featAcc(i) = 1;
    end
end
fprintf('\n');

ktr = gendat(tr, [5000 5000 5000]);
wrf = librandomforestc(ktr, 50);
tmp = +wrf;
featAcc(:, 5) = 0.5 - tmp.importance(:,4);

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

fprintf('Feature ranking and clustering:\n-----------------------------\n');
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

fprintf('Cluster Ranking:\n----------------\n');
for i = 1:nclust
    fprintf('Cluster %d, Accuracy %.3f\n', clustIdx(i), clustAcc(clustIdx(i)));
    fprintf('\t%s, ',clustFl{clustIdx(i)}{:});
    fprintf('\n');
end

%% What sort of acc do we get if we choose one feature from each of the 1st N clusters
load('G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\DataAllWin5NoBorder2.mat')
dataAll = changelablist(dataAll, 'Default')
dataAll = setprior(dataAll, 0);
global feats
feats = [9 15 20 23 7 6]; %ranked cluster
% feats = [16 23 21 25 38 20]; %wfso = featself(subData, libsvc([], proxm([], 'r', 2), 20), 6);
% feats = +wfso

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
% clfr.save('G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\configDTree.yaml');

clfrSvm = +w;
clfrSvm = +clfrSvm{2};
clfrSvm.save('G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\configSvm.yaml');

% save('G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\All_Svc_clfr3.mat', 'w', 'feats')
% save('G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\All_RandomForest_clfr_3class.mat', 'w', 'feats')
% save('G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\All_OpenCvRTrees_clfr_3class.mat', 'w', 'feats')
% save('G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\All_OpenCvDTree_clfr_3class.mat', 'w', 'feats')
% save('G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\All_OpenCvSvc_clfr_3class.mat', 'w', 'feats')
save('G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\All_OpenCvDTree_clfr_3class.mat', 'w', 'feats', 'subData')

s = [];
s.tr = single(+subData(:, feats));
s.trLab = single(+getnlab(subData));

%copy prev file contents before exec this line
cv.FileStorage('G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\trDataDTree.yaml', s);


%% How do other FS algorithms work in terms of what clusters they choose from
subData = gendat(dataAll, [2000 2000 2000]);
scaleData = subData*scalem(subData, 'variance');
[tr ts] = gendat(scaleData, 0.5);

[wfs rfs] = featself(subData, libsvc([], proxm([], 'r', 2), 20), 10);
% feats = +wfso
fidx = +wfso; % [16    23    21    25    38    20] - good 

% fl()
for i = 1:length(fidx)
    fprintf('Feature %d, %s, Cluster Id %d, Cluster Rank %d, Cluster Feats %s, \n', fidx(i), fl{fidx(i)}, lab(fidx(i)), find(clustIdx==lab(fidx(i))), [clustFl{lab(fidx(i))}{:}]);    
end

% wfso = featself(tr, naivebc, 6, ts);
%--------------------------------------------------------------------------
%% Compare 3 class, vs 2 class vs OCC classifiers
close all hidden; clear all;

load('G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\DataAllWin5NoBorder2.mat');
dataAll = changelablist(dataAll, 'Default');
dataAll = setprior(dataAll, 0);
subData = gendat(dataAll, [2000 4000 2000]);
subData = dataAll;

feats = [9 15 20 23 7 6]; %ranked cluster
% feats = [16 23 21 25 38 20]; %wfso = featself(subData, libsvc([], proxm([], 'r', 2), 20), 6);
ri = 1;
res(ri).data = subData;
res(ri).feats = feats;

%% 3 class
% [err, cerr, res(ri).nlabOut, tclassf, res(ri).tress] = prcrossval(res(ri).data(:, res(ri).feats), scalem([], 'variance')*libsvc([], proxm([], 'r', 2), 20), 5, 1);
ll = getlablist(subData);
[err, cerr, res(ri).nlabOut, tclassf, res(ri).tress] = prcrossval(res(ri).data(:, res(ri).feats), qdc*affine(diag([1 10 1]), [], ll, ll), 5, 1);

res(ri).c = confmat(getnlab(subData), res(ri).nlabOut);
res(ri).cn = res(ri).c./repmat(sum(res(ri).c, 2), 1, size(res(ri).c, 2))
1 - mean(diag(res(ri).cn)) %0.192788729769175
res(ri).crn = ReduceConfMat(res(ri).cn, {[1 3], [2]}, true)
1 - mean(diag(res(ri).crn)) %0.0948738817959547

out2Class(:,1) = max(res(ri).tress(:, [1 3]), [], 2);
out2Class(:,2) = res(ri).tress(:, 2);
otherIdx = getnlab(subData) == getclassi(subData, 'Background') | getnlab(subData) == getclassi(subData, 'Tree');
ll = cellstr(getlab(subData));
ll(otherIdx) = {'Other'};
out2Class = prdataset(out2Class, ll);
out2Class = setfeatlab(out2Class, {'Other', 'Spekboom'});
res(ri).out2Class = setprior(out2Class, 0);

res(ri).e = roc(res(ri).out2Class, 2, 100);

figure
plote(res(ri).e, 'b--')

subData_ = gendat(dataAll, [2000 3000 2000]);
subData_ = dataAll;
w = qdc*affine(diag([1 10 1]), [], ll, ll);
w = subData_(:, feats)*w;
subData_(:, feats)*w*testc
save('G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\All_Qdc_clfr_3class.mat', 'w', 'feats')


%% 2 class
ri = 2;
nlab = getnlab(subData);
otherIdx = getnlab(subData) == getclassi(subData, 'Background') | getnlab(subData) == getclassi(subData, 'Tree');
ll = cellstr(getlab(subData));
ll(otherIdx) = {'Other'};
subData = addlabels(subData, char(ll), '2Class');
subData = changelablist(subData, '2Class');
subData = setprior(subData, 0);
res(ri).data = subData;
res(ri).feats = feats;

[err, cerr, res(ri).nlabOut, tclassf, res(ri).tress] = prcrossval(res(ri).data(:, res(ri).feats), scalem([], 'variance')*libsvc([], proxm([], 'r', 2), 20), 5, 1);
% [err, cerr, res(ri).nlabOut, tclassf, res(ri).tress] = prcrossval(res(ri).data(:, res(ri).feats), qdc, 5, 1);
res(ri).c = confmat(getnlab(res(ri).data), res(ri).nlabOut);

res(ri).cn = res(ri).c./repmat(sum(res(ri).c, 2), 1, size(res(ri).c, 2))
1 - mean(diag(res(ri).cn)) %0.192788729769175

% [tr ts] = gendat(subDataAll, 0.5);
out = prdataset(res(ri).tress, getlab(res(ri).data));
out = setfeatlab(out, getlablist(res(ri).data));
res(ri).out = setprior(out, 0);
res(ri).e = roc(res(ri).out, 1, 100);

hold on;
plote(res(ri).e, 'r--')

%% occ
ri = 3;
subData = changelablist(subData, 'Default');
subDataOcc = oc_set(subData, 'Spekboom');
% subDataOcc = gendat(subDataOcc, [1000 1000]);
subDataOcc = setprior(subDataOcc, 0);
res(ri).data = subDataOcc;
% res(ri).feats = [16 23 21 25 38 20];
res(ri).feats = feats; %[16 23 21 25 38 20];

% [err, cerr, res(ri).nlabOut, tclassf, res(ri).tress] = prcrossval(res(ri).data(:, res(ri).feats), scalem([], 'variance')*svdd([], 0.1, 1), 5, 1);
% [err, cerr, res(ri).nlabOut, tclassf, res(ri).tress] = prcrossval(res(ri).data(:, res(ri).feats([1 2 3 5 6])), scalem([], 'variance')*knndd([], 0.1, 7), 5, 1); %good
% [err, cerr, res(ri).nlabOut, tclassf, res(ri).tress] = prcrossval(res(ri).data(:, res(ri).feats([1 2 3])), scalem([], 'variance')*newsvdd([], 0.1, 2), 5, 1);
% [err, cerr, res(ri).nlabOut, tclassf, res(ri).tress] = prcrossval(res(ri).data(:, res(ri).feats([1 2 3 5])), mcd_gauss_dd([], .1), 5, 1); %good
[err, cerr, res(ri).nlabOut, tclassf, res(ri).tress] = prcrossval(res(ri).data(:, res(ri).feats([1 2 3 4 5 6])), scalem([], 'variance')*mog_dd([], 0.1, [3 5], 'full'), 5, 1); %best
% [err, cerr, res(ri).nlabOut, tclassf, res(ri).tress] = prcrossval(res(ri).data(:, res(ri).feats([1 2 3 5])), scalem([], 'variance')*parzen_dd([], 0.1, 1.15), 5, 1); %okish

% [tr ts] = gendat(subDataOcc, 0.5);
% w = tr(:, feats([1 2 3 5]))*(scalem([], 'variance')*newsvdd([], 0.1, 2));
% w = tr(1:end, res(ri).feats([1 2 3 5]))*(scalem([], 'variance')*parzen_dd([], 0.1, 1.5));
% ts(:, res(ri).feats([1 2 3 5]))*w*testc

res(ri).c = confmat(getnlab(res(ri).data), res(ri).nlabOut);

res(ri).cn = res(ri).c./repmat(sum(res(ri).c, 2), 1, size(res(ri).c, 2))
1 - mean(diag(res(ri).cn)) 

% [tr ts] = gendat(subDataAll, 0.5);
out = prdataset(res(ri).tress, getlab(res(ri).data));
out = setfeatlab(out, getlablist(res(ri).data));
res(ri).out = setprior(out, 0);
res(ri).e = roc(res(ri).out, 2, 100);

hold on;
plote(res(ri).e, 'm--')

w = scalem([], 'variance')*mog_dd([], 0.1, [3 5], 'full')
w = subDataOcc(:, feats)*w
subDataOcc(:, feats)*w*testc

save('G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\All_Mogc_clfr_Occ.mat', 'w', 'feats')

%--------------------------------------------------------------------------
%% Make per area/habitat classifiers and test against all areas classifier
close all hidden; clear all;

load('G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\DataAllWin5NoBorder2.mat');
feats = [9 15 20 23 7 6]; %ranked cluster

dataAll = changelablist(dataAll, 'Area');
areas = {'Groenfontein', 'Matjiesvlei', 'Rooiberg', 'GrootKop'}; % cellstr(getlablist(dataAll));
areaLbls = cellstr(getlabels(dataAll));

res = [];
for i = 1:length(areas)
    fprintf('%s\n', areas{i});
    areaIdx = ~cellfun(@isempty, strfind(areaLbls, areas{i}));
    subData = dataAll(areaIdx,:);

    subData = changelablist(subData, 'Default');
    subData = setprior(subData, 0);
    subData = gendat(subData, [2000 4000 2000]);
    
    res(i).data = subData;
    [err, cerr, res(i).nlabOut, tclassf, res(i).tress] = prcrossval(res(i).data(:, feats), scalem([], 'variance')*librandomforestc([], 100), 5, 1);
%     [err, cerr, res(i).nlabOut, tclassf, res(i).tress] = prcrossval(res(i).data(:, feats), scalem([], 'variance')*libsvc([], proxm([], 'r', 2), 20), 5, 1);
% [err, cerr, res(ri).nlabOut, tclassf, res(ri).tress] = prcrossval(res(ri).data(:, res(ri).feats), qdc, 5, 1);

    res(i).c = confmat(getnlab(subData), res(i).nlabOut);
    res(i).cn = res(i).c./repmat(sum(res(i).c, 2), 1, size(res(i).c, 2))
    1 - mean(diag(res(i).cn)) %0.192788729769175
    res(i).crn = ReduceConfMat(res(i).cn, {[1 3], [2]}, true)
    1 - mean(diag(res(i).crn)) %0.0948738817959547

    clear out2Class;
    out2Class(:,1) = max(res(i).tress(:, [1 3]), [], 2);
    out2Class(:,2) = res(i).tress(:, 2);
    otherIdx = getnlab(subData) == getclassi(subData, 'Background') | getnlab(subData) == getclassi(subData, 'Tree');
    ll = cellstr(getlab(subData));
    ll(otherIdx) = {'Other'};
    out2Class = prdataset(out2Class, ll);
    out2Class = setfeatlab(out2Class, {'Other', 'Spekboom'});
    res(i).out2Class = setprior(out2Class, 0);

    res(i).e = roc(res(i).out2Class, 2, 100);

end

c = sum(cat(3, res(:).c), 3)
cn = res(i).c./repmat(sum(res(i).c, 2), 1, size(res(i).c, 2))
1 - mean(diag(res(i).cn)) %0.192788729769175
crn = ReduceConfMat(res(i).cn, {[1 3], [2]}, true)
1 - mean(diag(res(i).crn)) %0.0948738817959547

out2Class = [];
for i = 1:length(res)
    out2Class = [out2Class; res(i).out2Class];
end
e = roc(out2Class, 2, 100);

figure
plote(e)

for i = 1:length(res)
%     res(i).w = res(i).data(1:2:end, feats)*(scalem([], 'variance')*libsvc([], proxm([], 'r', 2), 20));
    res(i).w = res(i).data(1:2:end, feats)*(scalem([], 'variance')*librandomforestc([], 100));
    
    w = res(i).w;
    fileName = sprintf('G:/MSc GeoInformatics/Data/NGI/My Rectified/Ground Truth Images/%s_RandomForest_Clfr_3Class.mat', areas{i});
    res_ = res(i);
    save(fileName, 'w', 'feats', 'res_');
end

%--------------------------------------------------------------------------
%% Test clfr against jan vlok gt
% clear all;
 close all hidden; 
imageFileNames = {...
    'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321b_3172_12_0415_rgbn_XCALIB.tif';...
    'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321b_3172_12_0419_rgbn_XCALIB.tif';...
    'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321D_319_04_0121_rgbn_XCALIB.tif';...
    'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3322C_322_02_0056_rgbn_XCALIB.tif';...
    };
% outImageFileNames = {...
%     'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321b_3172_12_0415_rgbn_XCALIB_Out_AllQdc.tif';...
%     'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321b_3172_12_0419_rgbn_XCALIB_Out_AllQdc.tif';...
%     'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321D_319_04_0121_rgbn_XCALIB_Out_AllQdc.tif';...
%     'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3322C_322_02_0056_rgbn_XCALIB_Out_AllQdc.tif';...
%     };

gtFileNames = {...
    'D:\Data\Development\Projects\MSc GeoInformatics\Docs\Trimble\GroundTruthCombinedCorrected_TM21_WGS84.shp';...
    'D:\Data\Development\Projects\MSc GeoInformatics\Docs\Trimble\GroundTruthCombinedCorrected_TM21_WGS84.shp';...
    'D:\Data\Development\Projects\MSc GeoInformatics\Docs\Trimble\GroundTruthCombinedCorrected_TM21_WGS84.shp';...
    'D:\Data\Development\Projects\MSc GeoInformatics\Docs\Trimble\GroundTruthCombinedCorrected_TM23_WGS84.shp';...
    };

% clfrFileNames = {...
%     'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\All_Knnc_clfr.mat';...
%     'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\All_Knnc_clfr.mat';...
%     'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\All_Knnc_clfr.mat';...
%     'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\All_Knnc_clfr.mat';...
%     };
% global w feats
clfrFileNames = {... %empty to use workspace vars
    '';...
    '';...
    '';...
    '';...
    };

%% Extract features and classify for gt regions
% shapeFileName = 'D:\Data\Development\Projects\MSc GeoInformatics\Docs\Trimble\GroundTruthCombinedCorrected_TM21_WGS84.shp';
close all hidden;

resC = {};
res = [];
for g = 1:length(gtFileNames)
    resC{g} = ValidateClfrAgainstJvGroundTruth(imageFileNames{g}, gtFileNames{g}, clfrFileNames{g});
end

res = [resC{:}];

% res = res_;
%% Display Results
for i = 1:length(res)
      fprintf('%s %f %f\n', res(i).name, res(i).gtCover, res(i).clfCover);
end

fprintf('Canopy cover error Mean(abs) Std(abs): %f (%f)\n', mean(abs([res.gtCover]-[res.clfCover])), std(abs([res.gtCover]-[res.clfCover])));
fprintf('Canopy cover error Median Mad: %f (%f)\n', median([res.gtCover]-[res.clfCover]), mad([res.gtCover]-[res.clfCover]));

%% Save results
% fn = 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\All_RandomForest_clfr_3class.mat';
load(clfrFileNames{1})
save(clfrFileNames{1}, 'w', 'feats', 'res')

%% Make gt sub ims for Jan
close all hidden;

resI = {};
res = [];
for g = 1:length(gtFileNames)
%     resC{g} = ValidateClfrAgainstJvGroundTruth(imageFileNames{g}, gtFileNames{g}, clfrFileNames{g});
    resI{g} = MakeGroundTruthImage(imageFileNames{g}, gtFileNames{g}, clfrFileNames{g});
end

res = [resI{:}];

%%
close all hidden;
figure
f = 1;
for i = 1:length(res)
    subplot(2,2,f)
    imshow([res(i).rgbIm res(i).cirIm]);
    title(res(i).name);

    f = f + 1;
    if f > 4
        f = 1;
        figure;
    end
    fileIm = [res(i).rgbIm res(i).cirIm];
    fileName = ['D:\Data\Development\Projects\MSc GeoInformatics\Data\' res(i).name '.png'];
    imwrite(fileIm, fileName);
end

%--------------------------------------------------------------------------
%% Make output images

% close all; clear all;
warning off
% clfrFileNames = {...
%     'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\Groenfontein_RandomForest_Clfr_3Class.mat';...
%     'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\Matjiesvlei_RandomForest_Clfr_3Class.mat';...
%     'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\Rooiberg_RandomForest_Clfr_3Class.mat';...
%     'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\GrootKop_RandomForest_Clfr_3Class.mat';...
%     };
fn(1).dataFileName = 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321b_3172_12_0415_rgbn_data.mat';
fn(1).imageFileName = 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321b_3172_12_0415_rgbn_XCALIB.tif';
fn(1).outImageFileName = 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321b_3172_12_0415_rgbn_XCALIB_Out_DTree3.tif';
% fn(1).clfrFileName = 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321b_3172_12_0415_clfr.mat';
% fn(1).clfrFileName = 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\Groenfontein_RandomForest_Clfr_3Class.mat';
fn(1).clfrFileName = 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\All_Svc_clfr_3class.mat';

fn(2).dataFileName = 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321b_3172_12_0419_rgbn_data.mat';
fn(2).imageFileName = 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321b_3172_12_0419_rgbn_XCALIB.tif';
fn(2).outImageFileName = 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321b_3172_12_0419_rgbn_XCALIB_Out_DTree3.tif';
% fn(2).clfrFileName = 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321b_3172_12_0419_rgbn_clfr.mat';
% fn(2).clfrFileName = 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\Matjiesvlei_RandomForest_Clfr_3Class.mat';
fn(2).clfrFileName = 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\All_Svc_clfr_3class.mat';

fn(3).dataFileName = 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321D_319_04_0121_rgbn_data.mat';
fn(3).imageFileName = 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321D_319_04_0121_rgbn_XCALIB.tif';
fn(3).outImageFileName = 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321D_319_04_0121_rgbn_XCALIB_Out_DTree3.tif';
% fn(3).clfrFileName = 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321D_319_04_0121_clfr.mat';
% fn(3).clfrFileName = 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\Rooiberg_RandomForest_Clfr_3Class.mat';
fn(3).clfrFileName = 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\All_Svc_clfr_3class.mat';

fn(4).dataFileName = 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3322C_322_02_0056_rgbn_data.mat';
fn(4).imageFileName = 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3322C_322_02_0056_rgbn_XCALIB.tif';
fn(4).outImageFileName = 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3322C_322_02_0056_rgbn_XCALIB_Out_DTree3.tif';
% fn(4).clfrFileName = 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3322C_322_02_0056_clfr.mat';
% fn(4).clfrFileName = 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\GrootKop_RandomForest_Clfr_3Class.mat';
fn(4).clfrFileName = 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\All_Svc_clfr_3class.mat';

if (false)
% matlabpool close
for i = 1:length(fn)
    prload(fn(i).clfrFileName);
    blockproc(fn(i).imageFileName, [128 128], @(x)ClassifyIm(x, w, feats), 'Destination', ...
        fn(i).outImageFileName, 'UseParallel', true);    
end
else
    cd('D:\Data\Development\Projects\MSc GeoInformatics\Code\SpekboomClassifier\x64\Release');
    for i = 1:length(fn)
        fprintf('SpekboomClassifier.exe "%s" "%s" -o\n', fn(i).imageFileName, fn(i).outImageFileName);
    end
    
end

%%
for i = 1:length(fn)
    out = imread(fn(i).outImageFileName);

    [dum outC] = max(+out,[],3);
    outC = outC==2; %2;
    clear dum %out

    % outIm = reshape(+out(:,2), size(mask));
    im = imread(fn(i).imageFileName);
    %
    figure;
    h1 = subplot(1,3,1);
    imshow(out(:,:,2));
    h2 = subplot(1,3,2);
    imshow(im(:,:,[4 1 2])*16);
    h3 = subplot(1,3,3);
    imshow(im(:,:,[1 2 3])*16);
    linkaxes([h1 h2 h3], 'xy');    
end

%% Convert clfr output files to geotiffs
cd 'C:\OSGeo4W64\bin'

%fix no_data
for i = 1:length(fn)
    out = imread(fn(i).outImageFileName);
    gi = geotiffinfo(fn(i).imageFileName);
    newoutImageFileName = [fn(i).outImageFileName(1:end-4) '_GT.tif'];
    border = all(out==0, 3);
    out(out==0) = 1;
    out(repmat(border, [1 1 3])) = 0;
    geotiffwrite(newoutImageFileName, out, gi.SpatialRef, ...
        'GeoKeyDirectoryTag', gi.GeoTIFFTags.GeoKeyDirectoryTag);
    fprintf('.');
    %copy and paste into dos win
end
fprintf('\n');

%set spatial info
for i = 1:3
    newoutImageFileName = [fn(i).outImageFileName(1:end-4) '_GT.tif']; 
    dos(sprintf('gdal_edit.py -a_nodata 0 -a_srs "%s" "%s"\n', ... %
         'PROJCS["unnamed", GEOGCS["unknown", DATUM["unknown", SPHEROID["WGS84",6378137,298.257223563]], PRIMEM["Greenwich",0], UNIT["degree",0.0174532925199433]], PROJECTION["Transverse_Mercator"], PARAMETER["latitude_of_origin",0], PARAMETER["central_meridian",21], PARAMETER["scale_factor",1], PARAMETER["false_easting",0], PARAMETER["false_northing",0], UNIT["metre",1]]', ...
         newoutImageFileName));
end


for i = 4:4
    newoutImageFileName = [fn(i).outImageFileName(1:end-4) '_GT.tif'];
    dos(sprintf('gdal_edit.py -a_nodata 0 -a_srs "%s" "%s"\n', ... %-a_nodata 0 
         'PROJCS["unnamed", GEOGCS["unknown", DATUM["unknown", SPHEROID["WGS84",6378137,298.257223563]], PRIMEM["Greenwich",0], UNIT["degree",0.0174532925199433]], PROJECTION["Transverse_Mercator"], PARAMETER["latitude_of_origin",0], PARAMETER["central_meridian",23], PARAMETER["scale_factor",1], PARAMETER["false_easting",0], PARAMETER["false_northing",0], UNIT["metre",1]]', ...
         newoutImageFileName));
end

tmp = imread(newoutImageFileName);
imshow(tmp)

%------------------------------------------------------------------------
%% Test mexopencv and compare results with prtools

close all hidden; clear all;

load('G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\DataAllWin5NoBorder2.mat');
dataAll = changelablist(dataAll, 'Default');
dataAll = setprior(dataAll, 0);
subData = gendat(dataAll, [2000 4000 2000]);
% subData = dataAll;

feats = [9 15 20 23 7 6]; %ranked cluster

subData_ = subData*scalem(subData, 'variance'); %nb scale data
[tr ts] = gendat(subData_(:, feats), 0.5);

% [err, cerr, res(ri).nlabOut, tclassf, res(ri).tress] = prcrossval(res(ri).data(:, res(ri).feats), scalem([], 'variance')*libsvc([], proxm([], 'r', 2), 20), 5, 1);
% tr = testdatasize(tr,'objects');
% K = compute_kernel(tr,tr,proxm([], 'r', 2));
% K = min(K,K');   % make sure kernel is symmetric
% K = [[1:m]' K];  % as libsvm wants it
% 	                   % call libsvm

clfr = cv.SVM;

clfr.train(+tr, getnlab(tr), 'SVMType', 'C_SVC', 'KernelType', 'RBF', 'Gamma', 2, 'C', 15);

nlabOut = clfr.predict(+ts);

c = confmat(getnlab(ts), double(nlabOut));
cn = c./repmat(sum(c, 2), 1, size(c, 2))
1 - mean(diag(cn)) %0.192788729769175
crn = ReduceConfMat(cn, {[1 3], [2]}, true)
1 - mean(diag(crn)) %0.0948738817959547

w = scalem([], 'variance')*libsvc([], proxm([], 'r', 2), 15)
w = tr*w

ts*w*testc
c = confmat(getnlab(ts), ts*w*nlabeld);
cn = c./repmat(sum(c, 2), 1, size(c, 2))
1 - mean(diag(cn)) %0.192788729769175
crn = ReduceConfMat(cn, {[1 3], [2]}, true)
1 - mean(diag(crn)) %0.0948738817959547


clfr.save('G:/MSc GeoInformatics/Data/NGI/My Rectified/Ground Truth Images/config.yaml');

%%
% [err, cerr, nlabOut] = prcrossval(subData(:, feats), scalem([], 'variance')*librandomforestc([], 100), 5, 1);
[tr ts] = gendat(subData, 0.5);

c = cv.RTrees;
c.train(+tr, getnlab(tr), 'Priors', [1 1 1]/3, 'MaxNumOfTreesInTheForest', 100, 'NActiveVars', 7, 'CalcVarImportance', true, 'MaxDepth', 10, 'ForestAccuracy', 0.05);

nlabOut = c.predict(+ts);

c = confmat(getnlab(ts), double(nlabOut));
cn = c./repmat(sum(c, 2), 1, size(c, 2))
1 - mean(diag(cn))
crn = ReduceConfMat(cn, {[1 3], [2]}, true)
1 - mean(diag(crn))

w = librandomforestc([], 100);
w = tr*w

ts*w*testc
c = confmat(getnlab(ts), ts*w*nlabeld);
cn = c./repmat(sum(c, 2), 1, size(c, 2))
1 - mean(diag(cn)) %0.192788729769175
crn = ReduceConfMat(cn, {[1 3], [2]}, true)
1 - mean(diag(crn)) %0.0948738817959547

%% Test OpenCV feature extraction
close all hidden;
fileName = 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321b_3172_12_0415_rgbn_XCALIB.tif';

load('G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\DataAllWin5NoBorder2.mat');

load('G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\PcaAll.mat')
load('G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\PcaSpekBoom.mat')
load('G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\PcaRgG.mat')

feats = [9 15 20 23 7 6]; %ranked cluster
fl = cellstr(getfeatlab(dataAll));
fl(feats)

w = +wPcaAll;
pcaRgb = [w.rot;w.offset];
tmp = [dataAll(:,1:4),ones(size(dataAll,1),1)] * pcaRgb;
tmp2 = dataAll(:,1:4)*wPcaAll;

figure;
for i=1:4
    subplot(2,2,i)
    plot(+tmp(:,i))
    hold all
    plot(+tmp2(:,i))    
end

%% 
w = +wPcaRgG;
pcaRgg = [w.rot;w.offset];
rgG = dataAll(:,1:4)./repmat(sum(dataAll(:,1:4),2),1,4);
tmp = [rgG,ones(size(dataAll,1),1)] * pcaRgg;
tmp2 = dataAll(:,5:8)*wPcaRgG;

figure;
for i=1:4
    subplot(2,2,i)
    plot(+tmp(:,i))
    hold all
    plot(+tmp2(:,i))    
end

%%
clfr.save('G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\config.yaml');
s = [];
s.rgbPca = single(pcaRgb);
s.rggPca = single(pcaRgg);

cv.FileStorage('G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\config.yaml', s);

%%

pcaAll = cv.PCA(+dataAll(:, 1:4));
rgbPca = pcaAll.project(+dataAll(:, 1:4));

rgG = +dataAll(:, 1:4)./repmat(sum(+dataAll(:, 1:4), 2), 1, 4);
pcaRgG = cv.PCA(rgG);
rgGPca = pcaRgG.project(rgG);

% mean(mean(abs(rgbPca - +dataAll(:,15:18))))
pcaAll.eigenvalues
pcaAll.eigenvectors

[V,D] = eig(cov(+dataAll(:, 1:4)));

tmp = +wPcaAll
tmp.eigenvalues
tmp.rot
mean(+dataAll(:,1:4), 1)

featIm = ExtractFeaturesIm2(double(im), 'feats', feats);

im = imread(fileName);
im = im(1:4:end,1:4:end,:);

imr = reshape(im, [], size(im, 3));

imRgG = double(im)./(eps + repmat(sum(im,3), [1, 1, size(im, 3)]));
imRgGr = reshape(imRgG, [], size(im, 3));
pcaRgG = cv.PCA(imRgGr);
imPcaRgGr = pcaRgG.project(imRgGr);
imPcaRgG = reshape(imPcaRgGr, size(im));


figure;
for i = 1:size(imPca,3)
    subplot(2,2,i)
    imagesc(imPca(:,:,i))
    colormap gray;
end

figure;
for i = 1:size(imPcaRgG,3)
    subplot(2,2,i)
    imagesc(imPcaRgG(:,:,i))
    colormap gray;
end

%--------------------------------------------------------------------------
%% Make opencv config file containing clfr etc
close all hidden; clear all;

prload('G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\DataAllWin5NoBorder2.mat');

prload('G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\PcaAll.mat')
prload('G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\PcaSpekBoom.mat')
prload('G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\PcaRgG.mat')

% feats = [9 15 20 7 6]; %NB missing entropy!
feats = [9 15 20 23 7 6]; %ranked cluster
fl = cellstr(getfeatlab(dataAll));
fl(feats)

%make classifier
%subData = gendat(dataAll, )
dataAll = changelablist(dataAll, 'Default');
dataAll = setprior(dataAll, 0);
subData = gendat(dataAll(:, feats), [2000 4000 2000]);

wScale = scalem(subData, 'variance');
subDataScale = subData*wScale; %nb scale data

%%
% rfeats = [1:9 15:23];
% subData = gendat(dataAll, [27260 27260 27260])
[tr ts] = gendat(dataAll(:,feats), 0.5);
clfr = cv.NormalBayesClassifier;

clfr.train(+tr, getnlab(tr));
% % clfr.train(+tr, getnlab(tr), 'Priors', [1 1.5 1]/3, 'MaxNumOfTreesInTheForest', 10, 'NActiveVars', 4, 'CalcVarImportance', true, 'MaxDepth', 5, 'ForestAccuracy', 0.05);

nlabOut = clfr.predict(+ts);

c = confmat(getnlab(ts), double(nlabOut));
cn = c./repmat(sum(c, 2), 1, size(c, 2))
1 - mean(diag(cn))
crn = ReduceConfMat(cn, {[1 3], [2]}, true)
1 - mean(diag(crn))

clfr.save('G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\configNormalBayes.yaml');


%%
if false
    subData = gendat(dataAll(:,feats), [27260 27260 27260])
else
    subData = dataAll(:,feats);
end
[tr_ ts_] = gendat(subData, 0.5);
wScale = scalem(tr_, 'variance')
tr = tr_*wScale
ts = ts_*wScale

clfr = cv.SVM;

clfr.train((+tr), getnlab(tr), 'SVMType', 'C_SVC', 'KernelType', 'RBF', 'Gamma', 5, 'C', 1, 'ClassWeights', double([1; 1; 1.5]));

% clfr.train((+tr), getnlab(tr), 'SVMType', 'C_SVC', 'KernelType', 'RBF', 'Gamma', 5, 'C', 1, 'ClassWeights', double([1; 1; 1]));
% clfr.train((+tr), getnlab(tr), 'SVMType', 'C_SVC', 'KernelType', 'Sigmoid', 'Gamma', 5, 'coeff', 0, 'C', 1, 'ClassWeights', double([1; 1; 1]));

nlabOut = clfr.predict(+ts);

c = confmat(getnlab(ts), double(nlabOut));
cn = c./repmat(sum(c, 2), 1, size(c, 2))
1 - mean(diag(cn)) %0.192788729769175
crn = ReduceConfMat(cn, {[1 3], [2]}, true)
1 - mean(diag(crn)) %0.0948738817959547

% clfr.save('G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\configSvm.yaml');

%%
% rfeats = [1:9 15:23];
[tr ts] = gendat(dataAll(:,feats), 0.5);
clfr = cv.RTrees;
% clfr.train(+tr, getnlab(tr), 'Priors', [1 1 1]/3, 'MaxNumOfTreesInTheForest', 50, 'NActiveVars', 4, 'CalcVarImportance', false, 'MaxDepth', 7, 'ForestAccuracy', 0.05);
clfr.train(+tr, getnlab(tr), 'Priors', [1 1 1]/3, 'MaxNumOfTreesInTheForest', 10, 'NActiveVars', 4, 'CalcVarImportance', true, 'MaxDepth', 5, 'ForestAccuracy', 0.05);

nlabOut = clfr.predict(+ts);

c = confmat(getnlab(ts), double(nlabOut));
cn = c./repmat(sum(c, 2), 1, size(c, 2))
1 - mean(diag(cn))
crn = ReduceConfMat(cn, {[1 3], [2]}, true)
1 - mean(diag(crn))

% clfr.save('G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\configRtrees.yaml');

%%
% rfeats = [1:9 15:23];
[tr ts] = gendat(dataAll(:,feats), 0.5);
clfr = cv.DTree;
% clfr.train(+tr, getnlab(tr), 'Priors', [1 1 1]/3, 'MaxNumOfTreesInTheForest', 50, 'NActiveVars', 4, 'CalcVarImportance', false, 'MaxDepth', 7, 'ForestAccuracy', 0.05);
clfr.train(+tr, getnlab(tr), 'Priors', [1 1.5 1]/4.5, 'MaxDepth', 8, 'Use1seRule', false, ...
    'UseSurrogates', false, 'CVFolds', 5, 'TruncatePrunedTree', true, 'MinSampleCount', min(classsizes(tr))/20); %, 'CalcVarImportance', true, 'MaxDepth', 15, 'ForestAccuracy', 0.05);

nlabOut = clfr.predict(+ts);

c = confmat(getnlab(ts), double(nlabOut));
cn = c./repmat(sum(c, 2), 1, size(c, 2))
1 - mean(diag(cn))
crn = ReduceConfMat(cn, {[1 3], [2]}, true)
1 - mean(diag(crn))

clfr.save('G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\configDtree.yaml');


%%
[tr_ ts] = gendat(dataAll(:,feats), 0.5);
tr_ = gendat(tr_, [500 1000 500]);
wScale = scalem(tr_, 'variance')
tr = tr_*wScale
ts = ts*wScale

clfr = cv.KNearest;
clfr.train(+tr, getnlab(tr));
% clfr.train(+tr, getnlab(tr), 'Priors', [1 1.5 1]/3, 'MaxNumOfTreesInTheForest', 10, 'NActiveVars', 4, 'CalcVarImportance', true, 'MaxDepth', 5, 'ForestAccuracy', 0.05);

nlabOut = clfr.predict(+ts);

c = confmat(getnlab(ts), double(nlabOut));
cn = c./repmat(sum(c, 2), 1, size(c, 2))
1 - mean(diag(cn))
crn = ReduceConfMat(cn, {[1 3], [2]}, true)
1 - mean(diag(crn))

clfr.save('G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\configKNearest.yaml');

%%
w = +wPcaAll;
rgbPca = [w.rot;w.offset];
w = +wPcaRgG;
rggPca = [w.rot;w.offset];
w = +wScale;
scale = [w.rot;w.offset];

s = [];
s.rgbPca = single(rgbPca);
s.rggPca = single(rggPca);
s.scale = single(scale);
% s.tr = single(+tr_);
% s.trLab = single(+getnlab(tr));

%copy prev file contents before exec this line
cv.FileStorage('G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\config.yaml', s);

s = [];
s.tr = single(+tr_);
s.trLab = single(+getnlab(tr_));

%copy prev file contents before exec this line
cv.FileStorage('G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\trData.yaml', s);


%- wPcaAll2 is not the same as wPcaAll and so feats have changed

%---------------------------------------------------------------------------
%% Debug OpenCV results

fn = 'D:/results.yaml.gz';
s = cv.FileStorage(fn);

fid = fopen('D:/block.txt', 'r');
res = fscanf(fid,'%d,');
fclose(fid);

block = reshape(res, [4 128 128]);
block = shiftdim(block,1);

figure
imagesc(block(:,:,[4 1 2])/5107)

fid = fopen('D:/feats.txt', 'r');
resf = fscanf(fid,'%g,');
fclose(fid);

featImOpenCV = reshape(resf, [5 128 128]);
featImOpenCV = shiftdim(featImOpenCV,1);
PlotMultibandImage(featImOpenCV)
% feats = [9 15 20 23 7 6]; %ranked cluster

featImMatlab = ExtractFeaturesIm2(double(block)/(2^14), 'feats', [9 15 20 7 6]);
featImMatlab = featImMatlab(:, :, [9 15 20 7 6]);
PlotMultibandImage([featImOpenCV featImMatlab])
PlotMultibandImage([featImOpenCV - featImMatlab])

outOpenCv = clfr.predict(reshape(featImMatlab, [], 5));

figure;
imagesc(reshape(out,128, 128))

wSvc = libsvc([], proxm([], 'r', 5), 20);
wSvc = tr*wSvc;

f = prdataset(reshape(featImMatlab, [], 5))*m;
outMatlab = f(1:end/2,:)*wSvc*nlabeld;
outMatlab2 = f(end/2+1:end,:)*wSvc*nlabeld;
figure;
imagesc((reshape(+([outMatlab;outMatlab2]), 128, 128))')

figure;
imagesc([(s.entropy)' featImMatlab(:, :, 4)])

%-------------------------------------------------------------------------
%% Validation - get results for all jan vlok gt
%
% dataFileName = 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321b_3172_12_0415_rgbn_data.mat';
% imageFileName = 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321b_3172_12_0415_rgbn_XCALIB.tif';
% outImageFileName = 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321b_3172_12_0415_rgbn_XCALIB_Out.tif';
% clfrFileName = 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321b_3172_12_0415_clfr.mat';

% dataFileName = 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321b_3172_12_0419_rgbn_data.mat';
% imageFileName = 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321b_3172_12_0419_rgbn_XCALIB.tif';
% outImageFileName = 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321b_3172_12_0419_rgbn_XCALIB_Out.tif';
% clfrFileName = 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321b_3172_12_0419_rgbn_clfr.mat';

% dataFileName = 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321D_319_04_0121_rgbn_data.mat';
% imageFileName = 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321D_319_04_0121_rgbn_XCALIB.tif';
% outImageFileName = 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321D_319_04_0121_rgbn_XCALIB_Out.tif';
% clfrFileName = 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321D_319_04_0121_clfr.mat';

% dataFileName = 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3322C_322_02_0056_rgbn_data.mat';
% imageFileName = 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3322C_322_02_0056_rgbn_XCALIB.tif';
% outImageFileName = 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3322C_322_02_0056_rgbn_XCALIB_Out.tif';
% clfrFileNames = 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3322C_322_02_0056_clfr.mat';
close all; %clear all;
imageFileNames = {...
    'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321b_3172_12_0415_rgbn_XCALIB.tif';...
    'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321b_3172_12_0419_rgbn_XCALIB.tif';...
    'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321D_319_04_0121_rgbn_XCALIB.tif';...
    'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3322C_322_02_0056_rgbn_XCALIB.tif';...
    };
% outImageFileNames = {...
%     'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321b_3172_12_0415_rgbn_XCALIB_Out_DTree3.tif';...
%     'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321b_3172_12_0419_rgbn_XCALIB_Out_DTree3.tif';...
%     'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321D_319_04_0121_rgbn_XCALIB_Out_DTree3.tif';...
%     'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3322C_322_02_0056_rgbn_XCALIB_Out_DTree3.tif';...
%     };
outImageFileNames = {...
    'G:\MSc GeoInformatics\Data\NGI\Classification\3321b_3172_12_0415_rgbn_XCALIB_DTree_Out.tif';...
    'G:\MSc GeoInformatics\Data\NGI\Classification\3321b_3172_12_0419_rgbn_XCALIB_DTree_Out.tif';...
    'G:\MSc GeoInformatics\Data\NGI\Classification\3321D_319_04_0121_rgbn_XCALIB_DTree_Out.tif';...
    'G:\MSc GeoInformatics\Data\NGI\Classification\3322C_322_02_0056_rgbn_XCALIB_DTree_Out.tif';...
    };

gtFileNames = {...
    'D:\Data\Development\Projects\MSc GeoInformatics\Docs\Trimble\GroundTruthCombinedCorrected_TM21_WGS84.shp';...
    'D:\Data\Development\Projects\MSc GeoInformatics\Docs\Trimble\GroundTruthCombinedCorrected_TM21_WGS84.shp';...
    'D:\Data\Development\Projects\MSc GeoInformatics\Docs\Trimble\GroundTruthCombinedCorrected_TM21_WGS84.shp';...
    'D:\Data\Development\Projects\MSc GeoInformatics\Docs\Trimble\GroundTruthCombinedCorrected_TM23_WGS84.shp';...
    };

clfrFileNames = {...
    'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\All_RandomForest_clfr.mat';...
    'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\All_RandomForest_clfr.mat';...
    'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\All_RandomForest_clfr.mat';...
    'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\All_RandomForest_clfr.mat';...
    };

res = [];
count = 1;
for g = 1:length(gtFileNames)
    s = shaperead(gtFileNames{g});
    gi = geotiffinfo(imageFileNames{g});
%     im = imread(imageFileNames{g});
    out = imread(outImageFileNames{g});
%     [dum outC] = max(+out, [], 3);
    outC = out==255; %1;NB
    for i = 1:length(s)
        [row,col] = map2pix(gi.SpatialRef, [s(i).X], [s(i).Y]);
        if (min(row)>0 && max(row)<=size(out,1) && min(col)>0 && max(col)<=size(out,2))
            rowIdx = floor(min(row)):ceil(max(row));
            colIdx = floor(min(col)):ceil(max(col));

            idx = find(isnan(s(i).X));

    %         idx = 1:idx(1)-1;
            mask = false(size(out,1), size(out,2));
            startIdx = 1;
            for j = 1:length(idx)
                polyIdx = startIdx:idx(j)-1;
                mask = mask | poly2mask(col(polyIdx), row(polyIdx), size(out,1), size(out,2));
                startIdx = idx(j)+1;
            end
            fprintf('%s %d-%d, %f\n', s(i).Name, s(i).CoverMin, s(i).CoverMax, 100*sum(outC(mask))/sum(mask(:)));

            res(count).name = s(i).Name;
            res(count).gtCover = mean([s(i).CoverMin s(i).CoverMax]);
            res(count).clfCover = 100*sum(outC(mask))/sum(mask(:));
            count = count + 1;

            if false
                figure;
                h1 = subplot(2,2,1);
                imshow(mask(rowIdx, colIdx));
                h2 = subplot(2,2,2);
                imshow(outC(rowIdx, colIdx));
                h3 = subplot(2,2,3);
                imshow(im(rowIdx, colIdx, [1 2 3])*16);
                h4 = subplot(2,2,4);
                imshow(im(rowIdx, colIdx, [4 1 2])*16);
                linkaxes([h1 h2 h3 h4], 'xy')
            end
        end
    end
end

fprintf('Canopy cover error Mean(abs) Std(abs): %f (%f)\n', mean(abs([res.gtCover]-[res.clfCover])), std(abs([res.gtCover]-[res.clfCover])));
fprintf('Canopy cover error Median Mad: %f (%f)\n', median([res.gtCover]-[res.clfCover]), mad([res.gtCover]-[res.clfCover]));

%% Find accuracy of clfr against Jan's guess and Jan's guess against orig GT
xlsGtFile = 'D:\Data\Development\Projects\MSc GeoInformatics\Docs\My Docs\Field Trip 2012\SpekGuess.xls';
[num, txt, raw] = xlsread(xlsGtFile);
raw = raw(2:end,:);
for i = 1:length(res)
    idx = find(strcmpi(raw(:, 1), res(i).name));
    res(i).jvCover = mean([raw{idx,2:end}]);
    fprintf('%s gt: %2.1f, clf: %2.1f, jv: %2.1f \n', res(i).name, res(i).gtCover, res(i).clfCover, res(i).jvCover);
end

idx = 1:length(res);
idx = setdiff(1:length(res), strmatch('MatjiesVlei', {res.name}))

fprintf('Clf-JV canopy cover error Mean(abs) Std(abs): %f (%f)\n', mean(abs([res(idx).jvCover]-[res(idx).clfCover])), std(abs([res(idx).jvCover]-[res(idx).clfCover])));
fprintf('Clf-JV Canopy cover error Median Mad: %f (%f)\n', median([res(idx).jvCover]-[res(idx).clfCover]), mad([res(idx).jvCover]-[res(idx).clfCover]));
fprintf('GT-JV canopy cover error Mean(abs) Std(abs): %f (%f)\n', mean(abs([res(idx).gtCover]-[res(idx).jvCover])), std(abs([res(idx).gtCover]-[res(idx).jvCover])));
fprintf('GT-JV Canopy cover error Median Mad: %f (%f)\n', median([res(idx).gtCover]-[res(idx).jvCover]), mad([res(idx).gtCover]-[res(idx).jvCover]));
fprintf('Canopy cover error Mean(abs) Std(abs): %f (%f)\n', mean(abs([res(idx).gtCover]-[res(idx).clfCover])), std(abs([res(idx).gtCover]-[res(idx).clfCover])));
fprintf('Canopy cover error Median Mad: %f (%f)\n', median([res(idx).gtCover]-[res(idx).clfCover]), mad([res(idx).gtCover]-[res(idx).clfCover]));

cellArray = [{res.name}' {res.gtCover}' {res.clfCover}' {res.jvCover}'];
cellArray = [{'Area', 'Field Ground Truth', 'Classfier Estimate', 'JV Estimate'}; cellArray];
cellArray = [cellArray; {'Mean Abs Error (Std Dev)', '', sprintf('%.2f(%.2f)', mean(abs([res.gtCover]-[res.clfCover])), std(abs([res.gtCover]-[res.clfCover]))), ...
    sprintf('%.2f(%.2f)', mean(abs([res.gtCover]-[res.jvCover])), std(abs([res.gtCover]-[res.jvCover])))}];

%-------------------------------------------------------------------------
%% Compare OpenCV with prtools classifiers
close all; clear all;
prload('G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\DataAllWin5NoBorder2.mat');
dataAll = changelablist(dataAll, 'Default');
dataAll = setprior(dataAll, 0);
feats = [9 15 20 23 7 6]; %ranked cluster

%% SVC
subData = gendat(dataAll(:, feats), [2000 4000 2000]);

[err, cerr, nlabOut, tclassf, tress] = prcrossval(subData, scalem([], 'variance')*libsvc([], proxm([], 'r', 5), 10), 5, 1);

c = confmat(getnlab(subData), nlabOut);
cn = c./repmat(sum(c, 2), 1, size(c, 2))
1 - mean(diag(cn))

[CRn] = ReduceConfMat(cn, {[1 3], [2]}, true)
1 - mean(diag(CRn))

%
% subData = gendat(dataAll(:, feats), [2000 4000 2000]);
    
% [err, cerr, nlabOut, tclassf, tress] = prcrossval(subData, opencvrtreec([], {'Priors', [1 1 1]/3, 'MaxNumOfTreesInTheForest', 100}), 5, 1);
[err, cerr, nlabOut, tclassf, tress] = prcrossval(subData, scalem([], 'variance')*opencvsvc([], [], ...
    {'SVMType', 'C_SVC', 'KernelType', 'RBF', 'Gamma', 5, 'C', 10, 'ClassWeights', double([1; 1; 1])}), 5, 1);

c = confmat(getnlab(subData), nlabOut);
cn = c./repmat(sum(c, 2), 1, size(c, 2))
1 - mean(diag(cn))

[CRn] = ReduceConfMat(cn, {[1 3], [2]}, true)
1 - mean(diag(CRn))

%% Decision Tree
subData = gendat(dataAll(:, feats), [2000 4000 2000]);
    
[err, cerr, nlabOut, tclassf, tress] = prcrossval(subData, treec, 5, 1);

c = confmat(getnlab(subData), nlabOut);
cn = c./repmat(sum(c, 2), 1, size(c, 2))
1 - mean(diag(cn))

[CRn] = ReduceConfMat(cn, {[1 3], [2]}, true)
1 - mean(diag(CRn))

%
% subData = gendat(dataAll(:, feats), [2000 4000 2000]);
    
% [err, cerr, nlabOut, tclassf, tress] = prcrossval(subData, opencvrtreec([], {'Priors', [1 1 1]/3, 'MaxNumOfTreesInTheForest', 100}), 5, 1);
[err, cerr, nlabOut, tclassf, tress] = prcrossval(subData, opencvdtreec([], 12), 5, 1);

c = confmat(getnlab(subData), nlabOut);
cn = c./repmat(sum(c, 2), 1, size(c, 2))
1 - mean(diag(cn))

[CRn] = ReduceConfMat(cn, {[1 3], [2]}, true)
1 - mean(diag(CRn))

%% Random Forest
subData = gendat(dataAll(:, feats), [2000 4000 2000]);
    
[err, cerr, nlabOut, tclassf, tress] = prcrossval(subData, librandomforestc([], 100, 4), 5, 1);

c = confmat(getnlab(subData), nlabOut);
cn = c./repmat(sum(c, 2), 1, size(c, 2))
1 - mean(diag(cn))

[CRn] = ReduceConfMat(cn, {[1 3], [2]}, true)
1 - mean(diag(CRn))

%
% subData = gendat(dataAll(:, feats), [2000 4000 2000]);
    
% [err, cerr, nlabOut, tclassf, tress] = prcrossval(subData, opencvrtreec([], {'Priors', [1 1 1]/3, 'MaxNumOfTreesInTheForest', 100}), 5, 1);
[err, cerr, nlabOut, tclassf, tress] = prcrossval(subData, opencvrtreec([], [], {'Priors', [1 1 1]/3, ...
    'MaxNumOfTreesInTheForest', 100, 'NActiveVars', 4, 'CalcVarImportance', true, 'MaxDepth', 20, 'ForestAccuracy', 0.05}), 5, 1);

c = confmat(getnlab(subData), nlabOut);
cn = c./repmat(sum(c, 2), 1, size(c, 2))
1 - mean(diag(cn))

[CRn] = ReduceConfMat(cn, {[1 3], [2]}, true)
1 - mean(diag(CRn))




%%

lln = cellstr(getlablistnames(dataAll));

for i = 1:length(lln)
    dataAll = changelablist(dataAll, lln{i});
    lln{i}
    ll = getlablist(dataAll)
end

%% find ttl polygons & pixels of each class

clear all; close all;
shapeFileNames = {...
    'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321B_3172_12_0415_GroundTruth.shp';...
    'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321B_3172_12_0419_RGBN_GroundTruth.shp';...
    'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321B_319_04_0121_GroundTruth.shp';...
    'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3322C_322_02_0056_GroundTruth.shp';...
    };

classList = {'Spekboom', 'Background', 'Tree'};
classPolyCount = [0 0 0];
polyCount = 0;
for i = 1:length(shapeFileNames)
    s = shaperead(shapeFileNames{i});
    l = {s.Class};
    %classList
    for j = 1:length(classList)
        idx = strcmpi(classList{j}, l);
        classPolyCount(j) = classPolyCount(j) + sum(idx);
    end
    polyCount = polyCount + length(s);
end

