%% Feature Extraction
clear all; close all;

shapeFileNames = {...
    'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321B_3172_12_0415_GroundTruth.shp';...
    'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321B_3172_12_0419_RGBN_GroundTruth.shp';...
    'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321B_319_04_0121_GroundTruth.shp';...
    'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3322C_322_02_0056_GroundTruth.shp';...
    };

imageFileNames = {...
    'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321b_3172_12_0415_rgbn_XCALIB.tif';...
    'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321b_3172_12_0419_rgbn_XCALIB.tif';...
    'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321D_319_04_0121_rgbn_XCALIB.tif';...
    'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3322C_322_02_0056_rgbn_XCALIB.tif';...
    };

dataFileNames = {...
    'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321b_3172_12_0415_rgbn_rgbData.mat';
    'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321b_3172_12_0419_rgbn_rgbData.mat';
    'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321D_319_04_0121_rgbn_rgbData.mat';
    'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3322C_322_02_0056_rgbn_rgbData.mat';
    };

%% Extract RGBNIR features to get PCA mappings
dataAll = prdataset();
for i = 1:length(imageFileNames)
    fprintf('Processing %s\n', imageFileNames{i});
    s = shaperead(shapeFileNames{i});
    gi = geotiffinfo(imageFileNames{i});
    im = imread(imageFileNames{i});

    [data imData] = ExtractRgb(im, s, gi.SpatialRef);
    save(dataFileNames{i}, 'data')
    dataAll = [dataAll; data];
end

dataAll = changelablist(dataAll, 'Default');
scatterdui(dataAll);
save 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\DataAllRgb.mat' dataAll
%% Find PCA mappings
load('G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\DataAllWin5NoBorder2.mat');
dataAll = dataAll(:,1:4);

dataAll = setprior(dataAll, 0);
% dataAll = dataAll/(2^14); %scale 0-1 as in ExtractFeaturesIm2

wPcaAll = pcam(dataAll, 4);
scatterdui(dataAll*wPcaAll);

wAvgAll = klm(dataAll, 4);
scatterdui(dataAll*wAvgAll);

spekboomIdx = getclassi(dataAll, 'Spekboom');
spekboomData = seldat(dataAll, spekboomIdx);

wPcaSpekboom = pcam(spekboomData, 4);
scatterdui(dataAll*wPcaSpekboom);

dataRgG = dataAll./repmat(sum(+dataAll, 2), 1, 4);
wPcaRgG = pcam(dataRgG, 4);
scatterdui(dataRgG*wPcaRgG);

save 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\PcaAll2.mat' wPcaAll
save 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\PcaSpekBoom2.mat' wPcaSpekboom
save 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\PcaRgG2.mat' wPcaRgG

%% Extract all features
% dataAllRgb = dataAll;
prload 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\PcaAll.mat'
prload 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\PcaSpekBoom.mat'
prload 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\PcaRgG.mat'
dataAll = prdataset();
for i = 1:length(imageFileNames)
    fprintf('Processing %s\n', imageFileNames{i});
    s = shaperead(shapeFileNames{i});
    gi = geotiffinfo(imageFileNames{i});
    im = imread(imageFileNames{i});

    [data imData] = ExtractFeatures3(im, s, gi.SpatialRef, 'wPcaSpekboom', wPcaSpekboom, 'wPcaAll', wPcaAll, 'wPcaRgG', wPcaRgG, 'border', 2);
    dataC{i} = data;
    save(dataFileNames{i}, 'data')
    dataAll = [dataC{i}; dataAll];
end

%hack to concat data with empty slope labels
slopeLab = [];
dataAll_ = prdataset();
for i = 1:length(imageFileNames)
    dataC{i} = changelablist(dataC{i}, 'Slope');
    slopeLab = [slopeLab; getlab(dataC{i})];
    dataAll_ = [dataAll_; dellablist(dataC{i}, 'Slope')];
end
keepIdx = slopeLab=='N'|slopeLab=='S';
dataAll = addlabels(dataAll_(keepIdx, :), slopeLab(keepIdx), 'Slope');
dataAll = changelablist(dataAll, 'Slope');
% scatterd(dataAll(:, [2 4]), 'legend')
% scatterdui(dataAll);

save 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\DataAllWin5NoBorder2.mat' dataAll

%% Extract object features
% dataAllRgb = dataAll;
prload 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\PcaAll.mat'
prload 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\PcaSpekBoom.mat'
prload 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\PcaRgG.mat'
dataAll = prdataset();
for i = 1:length(imageFileNames)
    fprintf('Processing %s\n', imageFileNames{i});
    s = shaperead(shapeFileNames{i});
    gi = geotiffinfo(imageFileNames{i});
    im = imread(imageFileNames{i});
%     internal.images.isFigureAvailable = false;
    [data imData] = ExtractObjectFeatures(im, s, gi.SpatialRef, 'wPcaSpekboom', wPcaSpekboom, 'wPcaAll', wPcaAll, 'wPcaRgG', wPcaRgG);
%     save(dataFileNames{i}, 'data')
    dataAll = [dataAll; data];
end
%     internal.images.isFigureAvailable = true;

dataAll = changelablist(dataAll, 'Default');

save 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\DataAllObjectMean.mat' dataAll

scatterdui(dataAll);
fl = cellstr(getfeatlab(dataAll));
figure
subplot(2,2,1)
plotsig(dataAll(:, 51:56), 'legend')
title(fl{51})
subplot(2,2,2)
plotsig(dataAll(:, 57:62), 'legend')
title(fl{57})
subplot(2,2,3)
plotsig(dataAll(:, 63:68), 'legend')
title(fl{63})

%% Select object features and test clfr
clear all; close all;
load 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\DataAllObjectMean.mat'
dataAll = changelablist(dataAll, 'Default');
dataAll = setprior(dataAll, 0);

% n = max(classsizes(dataAll));
% rDataAll = prdataset();
% for i = 1:3
%     classData = seldat(dataAll, i);
%     classData = gendat(classData, n);
%     rDataAll = [rDataAll; classData];
% end
% 
% ci = getclassi2(rDataAll, 'Spekboom');
% dataAll = oc_set(rDataAll, ci);


fl = cellstr(getfeatlab(dataAll));

%exclude nan features
idx = isnan(+dataAll);
data(idx)=0;
% dataAll = dataAll(:, idx);

[tr ts] = gendat(dataAll, 0.5);

% [wo ro] = featselo(tr, naivebc, 4, ts);
% getfeatlab(ts*wo)
% feats = struct(wo).data{1};

% [wfs rfs] = featself(tr, scalem([], 'variance')*librandomforestc([], 100,1), 6, ts);
% [wfs rfs] = featself(tr, scalem([], 'variance')*libsvc([], proxm([], 'r',
% % 1)), 6, ts); feats =  [78    14     9    27    11 22];
% [wfs rfs] = featself(tr, qdc, 6, ts); %feats = [78     9    62     6    50    51];
[wfs rfs] = featself(tr, scalem([], 'variance')*knnc([], 3), 5, ts); %feats = [ 5     2    62    60    74     1];
getfeatlab(ts*wfs)
feats = struct(wfs).data{1};

% feats =  [78    14     9    27    11 22];

feats = [8     2    61    29    14];
% 
% [wbe rbe] = featselb(tr, qdc, 1, ts);
% getfeatlab(ts*wbe)
% feats = [35 61 80 40]; % 61 80 40];
% fl(feats)
% 
% feats = [6 9];% 51:56;
% w = tr(:, feats)*knnc([], 5);
w = tr(:, feats)*librandomforestc;

c = confmat(ts(:, feats)*w);
cn = c./repmat(sum(c, 2), 1, size(c, 2))
1 - mean(diag(cn))

[CRn] = ReduceConfMat(cn, {[1 3], [2]}, true)
1 - mean(diag(CRn))

[err, cerr, nlabOut] = prcrossval(dataAll(:, feats(1:5)), scalem([], 'variance')*knnc([], 3), 10, 1);
% [err, cerr, nlabOut] = prcrossval(dataAll(:, feats), qdc, 10, 1);
% [err, cerr, nlabOut] = prcrossval(dataAll, librandomforestc, 10, 1);
% [err, cerr, nlabOut] = prcrossval(dataAll(:, feats),  scalem([], 'variance')*libsvc([], proxm([], 'r', 1)), 10, 1);

c = confmat(getnlab(dataAll), nlabOut);
cn = c./repmat(sum(c, 2), 1, size(c, 2))
1 - mean(diag(cn))

[CRn] = ReduceConfMat(cn, {[1 3], [2]}, true)
1 - mean(diag(CRn))

% NOTES
% ----------
% - Feature selection almost always chooses texture or sliding win stats
% even though they don't look great in the scatter
% - This changes a lot with feat sel alg and criterion used
% - This suggests there is something in the texture even though it doesn't
% scatter nicely
% - For some (elusive) combinations of feature and clfr, the performance is
% very good and significantly (~2x) better than per-pixel.  I would be careful
% about reading too much into this though, there are many feats and few
% objects!
% - Bear in mind that crossval is important with such a small dataset
% otherwise results on a single test set give lots of variation.

%% Select per-pixel features and test clfr to compare to per-object 
close all hidden; clear all;
load 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\DataAllWin5NoBorder.mat';
dataAll = setprior(dataAll, 0);

% n = min(classsizes(dataAll));
% rDataAll = prdataset();
% for i = 1:3
%     classData = seldat(dataAll, i);
%     classData = gendat(classData, n);
%     rDataAll = [rDataAll; classData];
% end
% 
% ci = getclassi2(rDataAll, 'Spekboom');
% dataAll = oc_set(rDataAll, ci);

fl = cellstr(getfeatlab(dataAll));

%exclude nan features
idx = isnan(+dataAll);
dataAll(idx)=0;
% dataAll = dataAll(:, idx);

[tr ts] = gendat(dataAll, 0.5);

% [wo ro] = featselo(tr, naivebc, 4, ts);
% getfeatlab(ts*wo)
% feats = struct(wo).data{1};

[wfs rfs] = featself(tr, qdc, 5, ts); %feats = [31     4    29     6    13];
[wfs rfs] = featself(tr(1:10:end,:), scalem([], 'variance')*knnc([], 9), 5, ts(1:10:end,:)); %feats = [31    17    38     7    12];
[wfs rfs] = featself(tr(1:50:end,:), scalem([], 'variance')*libsvc([], proxm([], 'r', 1)), 5, ts(1:50:end,:)); %feats = [ 25    18    29    39     7];
getfeatlab(ts*wfs)
feats = struct(wfs).data{1};
% 
% prmemory(10000000)
wsel = tr(1:10:end, :)*librandomforestc;
importance = struct(wsel).data{1}.importance
[dum idx] = sort(importance(:,5), 'descend');
fl(idx)
out = ts*wsel
c = confmat(getnlab(ts), out*nlabeld);
cn = c./repmat(sum(c, 2), 1, size(c, 2))
1 - mean(diag(cn)) %0.192788729769175
[CRn] = ReduceConfMat(cn, {[1 3], [2]}, true)
1 - mean(diag(CRn)) %0.0948738817959547


[wfs rfs] = featself(tr(1:100:end, :), scalem([], 'variance')*libsvc([], proxm([], 'r', 1)), 5, ts(1:100:end, :));
getfeatlab(ts*wfs)
feats = struct(wfs).data{1};
feats = [24    18     6     7    37];  %scalem([], 'variance')*libsvc([], proxm([], 'r', 1))


[err, cerr, nlabOut] = prcrossval(dataAll(1:50:end, feats(1:5)), scalem([], 'variance')*libsvc([], proxm([], 'r', 1)), 10, 1); %0.0772505691478417
[err, cerr, nlabOut] = prcrossval(dataAll(1:50:end, :), librandomforestc, 10, 1);
c = confmat(getnlab(dataAll(1:50:end, feats(1:5))), nlabOut);


[err, cerr, nlabOut] = prcrossval(dataAll(:, feats), qdc, 10, 1);
c = confmat(getnlab(dataAll(:, feats)), nlabOut);

[err, cerr, nlabOut] = prcrossval(dataAll(1:10:end, feats), scalem([], 'variance')*knnc([], 9), 10, 1);
c = confmat(getnlab(dataAll(1:10:end, feats)), nlabOut);

cn = c./repmat(sum(c, 2), 1, size(c, 2))
1 - mean(diag(cn)) %0.192788729769175

[CRn] = ReduceConfMat(cn, {[1 3], [2]}, true)
1 - mean(diag(CRn)) %0.0948738817959547

%% train and save some clfr's

subDataAll = gendat(dataAll, [9999 9999 9999]);
[tr ts] = gendat(subDataAll, 0.5);
[wfs rfs] = featself(tr(1:10:end,:), scalem([], 'variance')*knnc([], 9), 5, ts(1:10:end,:)); %feats = [ 25     7    38    11    29];
[wfs rfs] = featself(tr(1:10:end,:), scalem([], 'variance')*libsvc([], proxm([], 'r', 1)), 5, ts(1:10:end,:)); %feats = [ 25     2    23     7     6];
getfeatlab(ts*wfs)
feats = struct(wfs).data{1};

[err, cerr, nlabOut] = prcrossval(subDataAll(:, feats),  scalem([], 'variance')*knnc([], 9), 10, 1);
[err, cerr, nlabOut] = prcrossval(subDataAll(1:5:end, feats),  scalem([], 'variance')*libsvc([], proxm([], 'r', 1)), 10, 1);
[err, cerr, nlabOut] = prcrossval(subDataAll(1:5:end, :),  librandomforestc, 10, 1);
c = confmat(getnlab(subDataAll(1:5:end, :)), nlabOut);

cn = c./repmat(sum(c, 2), 1, size(c, 2))
1 - mean(diag(cn)) %0.192788729769175

[CRn] = ReduceConfMat(cn, {[1 3], [2]}, true)
1 - mean(diag(CRn)) %0.0948738817959547


subSubDataAll = subDataAll(1:10:end, feats);
w = scalem([], 'variance')*knnc([], 9);
% J = edicon(distm(subSubDataAll));
% subSubDataAll = subSubDataAll(J,:);

w = subSubDataAll*w;
save('G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\All_Knnc_clfr.mat', 'w', 'feats')

w = scalem([], 'variance')*libsvc([], proxm([], 'r', 1))
w = subDataAll(1:5:end,feats)*w
save('G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\All_Svc_clfr.mat', 'w', 'feats')

% subSubDataAll = subDataAll(1:10:end, :);
w = librandomforestc
w = subDataAll(1:5:end, :)*w
feats = 1:size(subDataAll, 2);
save('G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\All_RandomForest_clfr.mat', 'w', 'feats')

%% libsvc experimentation

close all hidden; clear all;
prwarning off

load 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\DataAllWin5NoBorder.mat';
dataAll = setprior(dataAll, 0);
fl = cellstr(getfeatlab(dataAll));
dataAll(isnan(+dataAll))=0;

subDataAll = gendat(dataAll, [1000 1000 1000]);
[tr ts] = gendat(subDataAll, 0.5);

help featselm; help feateval

% [wFsNn rFsNn] = featselo(tr, 'NN', 5, ts);
% fl(+wFsNn)

ws = tr*scalem([], 'variance')
[wFsKnn rFsKnn] = featself(tr*ws, 'NN', 5, ts*ws); %feats = [21    26    33    40    16];%feats = [ 25     2    23     7     6];
feats = [ 25     2    23     7     6];
fl(+wFsKnn)
% [wFsMaha rFsMaha] = featselo(tr, 'maha-m', 5, ts); %SLOW
% fl(+wFsMaha)
[wFsQdc rFsQdc] = featself(tr, qdc, 5, ts);
fl(+wFsQdc)

[wFsSvc rFsSvc] = featself(tr, scalem([], 'variance')*libsvc([], proxm([], 'r', 2), 1), 5, ts); %feats = [ 25    18    29    39     7]; prev val
fl(+wFsSvc)
feats = +wFsSvc 

%%
%get kernel param
subDataAll = gendat(dataAll, [2000 2000 2000]);
%priors are balanced, separate into 2 classes before first training
% otherIdx = getnlab(subDataAll)==getclassi(dataAll, 'Background') | getnlab(subDataAll)==getclassi(dataAll, 'Tree');
% ll = cellstr(getlab(subDataAll));
% ll(otherIdx) = {'Other'};
% subDataAll = addlabels(subDataAll, char(ll), '2Class');
% subDataAll = changelablist(subDataAll, '2Class');
% subDataAll = setprior(subDataAll, 0);

% [tr ts] = gendat(gendat(subDataAll, [100 100]), 0.5);

[wRbLibSvc, kernel, nu] = rblibsvc(tr(:, feats)*scalem(tr(:, feats), 'variance')); 
+kernel %2-4

%get C param (NB this changes with the size of the training data!!!)
reg = nan;
scaledData = subDataAll(:, feats)*scalem(subDataAll(:, feats), 'variance');
w = regoptc(scaledData, 'libsvc' ,{proxm([], 'r', 2), reg}, {proxm([], 'r', 2), 1}, [2], [0 0; 1 200], testc([],'soft'));
+w %look at '-c ##'  ~20-100

[wFsSvc rFsSvc] = featself(tr, scalem([], 'variance')*libsvc([], proxm([], 'r', 2), 20), 5, ts); %feats = [ 25    18    29    39     7]; prev val
fl(+wFsSvc)
feats = +wFsSvc 

[err, cerr, nlabOut, tclassf, tress] = prcrossval(subDataAll(:, feats),  scalem([], 'variance')*libsvc([], proxm([], 'p', 3), 20), 5, 1);
c = confmat(getnlab(subDataAll), nlabOut);

cn = c./repmat(sum(c, 2), 1, size(c, 2))
1 - mean(diag(cn)) %0.192788729769175

[CRn] = ReduceConfMat(cn, {[1 3], [2]}, true)
1 - mean(diag(CRn)) %0.0948738817959547

% [tr ts] = gendat(subDataAll, 0.5);
out = prdataset(tress, getlab(subDataAll));
out = setfeatlab(out, getlablist(subDataAll));
e = roc(out, 3, 100);

hold on;
plote(e, 'b--')

%% Look at signatures (rgG, LBP, pca)
figure;
subplot(2,3,1)
plotsig(subDataAll(1:10:end, 1:4), 'legend')
title('RGBIR')
subplot(2,3,2)
plotsig(subDataAll(1:10:end, 5:8), 'legend')
title('rgGIR')
subplot(2,3,3)
plotsig(subDataAll(1:10:end, 11:14), 'legend')
title('TC')
subplot(2,3,4)
plotsig(subDataAll(1:10:end, 15:18), 'legend')
title('PC')
subplot(2,3,5)
plotsig(subDataAll(1:10:end, 19:22), 'legend')
title('RC')



% ------------------------------------------------------------------------
%% CART

close all hidden; clear all;
prwarning off

load 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\DataAllWin5NoBorder.mat';
dataAll = setprior(dataAll, 0);
fl = cellstr(getfeatlab(dataAll));
dataAll(isnan(+dataAll))=0;

subDataAll = gendat(dataAll, [1000 1000 1000]);
[tr ts] = gendat(subDataAll, 0.5);

w = statsdtc(ts);

statsTree = getdata(w);

fi = statsTree.predictorImportance();
[dum fIdx] = sort(fi, 'descend');

fl(fIdx)
%pc4 & tc4 - seriously???

%%
% set(0,'RecursionLimit',1000)
[wfs rfs] = featself(tr, dtc, 5, ts);
fl(+wfs)


%% Random forest vs DTC and Random Forest regoptc
wrf = librandomforestc(tr)
ts*wrf*testc

wdtc = dtc(tr)
ts*wdtc*testc

% scaledData = subDataAll(:, feats)*scalem(subDataAll(:, feats), 'variance');
ntree = nan;
w = regoptc(subDataAll, 'librandomforestc' , {ntree}, {100}, [1], [50 60], testc([],'soft'));
+w %ntree = 111; NB dependent on number of objects

ntree = 100;
mtry = nan;
w = regoptc(subDataAll, 'librandomforestc' , {ntree, mtry}, {100, 5}, [2], [100 100;5 20], testc([],'soft'));
+w %mtry = 14; NB dependent on number of objects


%% compare the above
w = tr*librandomforestc;
ts*w*testc

w = tr*librandomforestc([], 100, 14);
ts*w*testc
%makes bugger all difference

% ------------------------------------------------------------------------
%% Visualise per-pixel features
clear all;close all;
load('G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\DataAllWin5NoBorder2.mat');
data = dataAll;
fl = cellstr(getfeatlab(data));
feats = [6 9];
% feats = [5 6];
fl(feats)

data = changelablist(data, 'Default');
classes = cellstr(getlablist(data));
for i = 1:length(classes)
    classIdx{i} = (getnlab(data)==getclassi(data, classes{i}));
end

data = changelablist(data, 'Area');
areas = cellstr(getlablist(data));

figure;
cols = {'r','g','b','c','m','y','k'};
icons = {'x','o','<'};
leg = {};
h = [];
% clear areaIdx;
for i = 1:length(areas)
    classIdx_ = getclassi(data, areas{i});
    areaIdx{i} = (getnlab(data) == classIdx_(1)); %beware (1) - is hack to deal with substring matching issues
    for j = 1:length(classes)
        idx = find(areaIdx{i} & classIdx{j});
        if (~isempty(idx))
            h_ = plot(+data(idx, feats(1)), +data(idx, feats(2)), [cols{rem(i, length(cols))} icons{j}]);
            h(end+1) = h_(1)
            leg{end+1} = [areas{i} '-' classes{j}];
            hold on
        end
    end
end
legend(h, leg);

figure;
for i = 1:length(areas)
    data = changelablist(data, 'Area');   
    if (i == 1)
        areaIdx{i} = (getnlab(data) == 1);
    else
        areaIdx{i} = (getnlab(data) == getclassi(data, areas{i}));
    end
    data = changelablist(data, 'Default');    
    a(i) = subplot(3,3,i);
    if (i == 1)
        scatterd2(data(areaIdx{i}, feats), 'legend')
    else
        scatterd2(data(areaIdx{i}, feats))
    end
    title(areas{i})
    grid on
end
linkaxes(a, 'xy');


data = changelablist(data, 'Default');    
figure
scatterdui(data)

%% Compare per-pixel window sizes
close all; clear all;
dataFileNames = {...
    'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\DataAllWin3.mat';...
    'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\DataAllWin5.mat';...
    'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\DataAllWin7.mat';...
    'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\DataAllWin9.mat';...
    'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\DataAllWin11.mat';...
    'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\DataAllWin13.mat';...
    };
load(dataFileNames{3});
dataAll = changelablist(dataAll, 'Default');
dataAll = setprior(dataAll, 0);
dataAll = dataAll(:,15:end);  %remove per-pixel features

n = min(classsizes(dataAll));
rDataAll = prdataset();
for i = 1:3
    classData = seldat(dataAll, i);
    classData = gendat(classData, n);
    rDataAll = [rDataAll; classData];
end

ci = getclassi2(rDataAll, 'Spekboom');
dataAll = oc_set(rDataAll, ci);
[tr ts] = gendat(dataAll, 0.5);
fl = cellstr(getfeatlab(dataAll));

%Find best window features

% [wfs rfs] = featselo(tr(1:10:end,:), qdc, 4,  ts(1:10:end,:));
% [wfs rfs] = featselo(tr(1:10:end,:), qdc, 4,  ts(1:10:end,:));
[wfs rfs] = featselo(tr(1:5:end,:), gauss_dd([], .1), 2,  ts(1:5:end,:));
getfeatlab(ts*wfs)

%%
feats = [17 4 15 24];
feats = [9 25 4 18];
feats = [9 25 4 18];
feats = [17];

for i = 1:length(dataFileNames)
    fprintf('Processing %s\n', dataFileNames{i});
    load(dataFileNames{i});
    dataAll = changelablist(dataAll, 'Default');
    dataAll = setprior(dataAll, 0);    
    dataAll = dataAll(:,15:end);  %remove per-pixel features    
    
    n = min(classsizes(dataAll));
    rDataAll = prdataset();
    for j = 1:3
        classData = seldat(dataAll, j);
        classData = gendat(classData, n);
        rDataAll = [rDataAll; classData];
    end
    
    ci = getclassi2(rDataAll, 'Spekboom');
    dataAll = oc_set(rDataAll, ci);
    [tr ts] = gendat(dataAll, 0.5);
    
    w = tr(1:10:end, feats)*gauss_dd([], .2);
    c = confmat(ts(1:10:end, feats)*w);
    cn = c./repmat(sum(c, 2), 1, size(c, 2));
%     [CRn] = ReduceConfMat(cn, {[1 3], [2]}, true);
    err(i,1) = 1 - mean(diag(cn));
%     err(i,2) = 1 - mean(diag(CRn));    
end

figure;
plot([3 5 7 9 11 13], err,'-o');
grid on;
legend({'3 class', '2 class'})

% NOTES
%-------------------------------------------------------------------------
%- None of the texture features (except perhaps entropyIrRat) appear very useful
%- Some of the sliding win features (mean, median) are good
%- The accuracy vs win size characteristic varies a lot with features 
%and classifier used.  5 is probably best.

%% Make feature images for inspection in ArcMap
close all; clear all;
imageFileNames = {...
    'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321b_3172_12_0415_rgbn_XCALIB.tif';...
    'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321b_3172_12_0419_rgbn_XCALIB.tif';...
    'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321D_319_04_0121_rgbn_XCALIB.tif';...
    'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3322C_322_02_0056_rgbn_XCALIB.tif';...
    };
% matlabpool
cd 'C:\OSGeo4W64\bin'
numBands = 9;
outBandFileName = {};
for i = 1:length(imageFileNames)
    fprintf('Processing %s\n', imageFileNames{i});
    in = imfinfo(imageFileNames{i});
    outFileNames{i} = [imageFileNames{i}(1:end-4) '_FEATS.tif'];
    writeAdapter = MultiBandTiffAdapter(outFileNames{i}, 'w', [in.Height in.Width numBands] , double(0));
%     blockproc(imageFileNames{i}, [256 256], @(x)ExtractFeaturesIm2(x), 'Destination', ...
%         [imageFileNames{i}(1:end-4) '_FEATS.tif'], 'UseParallel', false);
    blockproc(imageFileNames{i}, [512 512], @(x)ExtractFeaturesIm2(x), 'Destination', ...
       writeAdapter, 'UseParallel', false);
    fprintf('\n');
    writeAdapter.close();
%    ts = Tiff(imageFileNames{i});
%    td = Tiff(outFileNames{i}, 'r+');   
    for j = 1:numBands
        outBandFileName{i, j} = sprintf('%s_Band%d.tif', outFileNames{i}(1:end-4), j);
        dos(sprintf('gdalcopyproj.py "%s" "%s"\n', [outFileNames{i}(1:end-10) '.tif'], outBandFileName{i, j}));
        dos(sprintf('gdal_edit.py -a_nodata 0 "%s"\n', outBandFileName{i, j}));
    end
end

for i = 1:length(outFileNames)
    bandFileNames = sprintf('"%s" ', outBandFileName{i,:});
    bandFileNames(end)=[];
    dos(sprintf('gdal_merge.py -of GTiff -co "BIGTIFF=YES" -a_nodata 0 -o "%s" -separate %s', outFileNames{i}, bandFileNames));
%     dos(sprintf('gdal_merge.py -o "%s" %s', outFileNames{i}, bandFileNames));
end

%% Feature selection
% 

% load(dataFileName);
data = dataAll;
data = changelablist(data, 'Default')
% p = [1 50 1];
% data = setprior(data, p./sum(p));
data = setprior(data, 0);
% [tr ts] = gendat(data(1:10:end,:), 0.5);
[tr ts] = gendat(data, 0.5);
fl = cellstr(getfeatlab(data));


% [wbe, r] = featselb(tr(1:10:end,:), scalem([], 'variance')*naivebc, 1, ts(1:10:end,:));
[wbe, rbe] = featselb(tr(1:10:end,:), qdc, 1, ts(1:10:end,:));

figure;
h = plot(42:-1:1, rbe(:,2), 'ko-');
fl(abs(rbe(2:end,3)))

[wfs, rfs] = featself(tr(1:10:end,:), qdc, 20, ts(1:10:end,:));

figure;
h = plot(rfs(:,1), rfs(:,2), 'ko-');
fl(abs(rfs(:,3)))

[wo ro] = featselo(tr(1:10:end,:), qdc, 4,  ts(1:10:end,:));
getfeatlab(ts*wo)

% set(gca, 'XTick', 1:size(r, 1));
% set(gca, 'XTickLabel', 1:size(dataAll, 2));

% scatterdui(tr)

% [wfs,r] = featselb(tr, naivebc, 3, ts);
% getfeatlab(tr*wfs) 
%---56
% R           
% G           
% stdI
%
% rN          
% bN          
% entropyIrRat
%---56

%---121
% nirN        
% G           
% bN    
%
% nirN        
% G           
% bN    
%---121

% featselo
% G           
% R           
% nirN 

% featself
% NDVI        
% R           
% entropyIrRat

% featselb
% bN          
% irRat       
% entropyIrRat

% featself(parzenc*scalem([], 'variance'))
% irRat       
% bN          
% NDVI      

% featself(tr, pca([], .9)*qdc, 3, ts);
% NDVI        
% rN          
% bN 

% featself(naivebc)
% irRat       
% bN          
% NDVI        

feats = [1 2];
% feats = [2 7 8];
% feats = [4 9];
% feats = [6 9];
% feats = [7 10];
% feats = [2 3 10];
w = qdc(tr(:, feats)) %*classReduceM (a, classIdxCell)

c = confmat(ts(:, feats)*w);
cn = c./repmat(sum(c, 2), 1, size(c, 2))
[CRn] = ReduceConfMat(cn, {[1 3], [2]}, true)
1 - mean(diag(CRn))

save(clfrFileName, 'w', 'feats');

figure;
scatterd2(tr(:, feats), 'legend')
plotc(w)


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
close all; clear all;
imageFileNames = {...
    'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321b_3172_12_0415_rgbn_XCALIB.tif';...
    'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321b_3172_12_0419_rgbn_XCALIB.tif';...
    'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321D_319_04_0121_rgbn_XCALIB.tif';...
    'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3322C_322_02_0056_rgbn_XCALIB.tif';...
    };
outImageFileNames = {...
    'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321b_3172_12_0415_rgbn_XCALIB_Out_AllQdc.tif';...
    'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321b_3172_12_0419_rgbn_XCALIB_Out_AllQdc.tif';...
    'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321D_319_04_0121_rgbn_XCALIB_Out_AllQdc.tif';...
    'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3322C_322_02_0056_rgbn_XCALIB_Out_AllQdc.tif';...
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
%% Extract features and classify for gt regions
% shapeFileName = 'D:\Data\Development\Projects\MSc GeoInformatics\Docs\Trimble\GroundTruthCombinedCorrected_TM21_WGS84.shp';
close all hidden;

res = [];
for g = 1:length(gtFileNames)
    res = [res, ValidateClfrAgainstJvGroundTruth(imageFileNames{g}, gtFileNames{g}, clfrFileNames{g})];
end

for i = 1:length(res)
      fprintf('%s %f %f\n', res(i).name, res(i).gtCover, res(i).clfCover);
%     for j = 1:length(res(i).name)
%         resC(count)
%         count = count+1;
%         fprintf('%s %f %f\n', res(i).name{j}, res(i).gtCover(j), res(i).clfCover(j));
%     end
end

fprintf('Canopy cover error: %f (%f)\n', mean(abs([res.gtCover]-[res.clfCover])), std(abs([res.gtCover]-[res.clfCover])));

%% Extract classification results from previously produced clfr output images
count = 1;
for g = 1:length(gtFileNames)
    s = shaperead(gtFileNames{g});
    gi = geotiffinfo(imageFileNames{g});
%     im = imread(imageFileNames{g});
    out = imread(outImageFileNames{g});
    [dum outC] = max(+out, [], 3);
    outC = outC==2; %1;NB
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

            name{count} = s(i).Name;
            gtCover(count) = mean([s(i).CoverMin s(i).CoverMax]);
            clfCover(count) = 100*sum(outC(mask))/sum(mask(:));
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

fprintf('Canopy cover error: %f (%f)\n', mean(abs(gtCover-clfCover)), std(abs(gtCover-clfCover)));
name = name';
res = name;
res = [res num2cell(gtCover(:)) num2cell(clfCover(:))];

%%
%convert output files to geotiffs
cd 'C:\OSGeo4W64\bin'
% gi.SpatialRef
% geotiffwrite([outImageFileNames{1}(1:end-4) '_GT.tif'], out, gi.SpatialRef, ...
%     'GeoKeyDirectoryTag', gi.GeoTIFFTags.GeoKeyDirectoryTag);
% dos(sprintf('gdal_edit.py -a_nodata 0 -a_srs "%s" "%s"', ...
%     'PROJCS["unnamed", GEOGCS["unknown", DATUM["unknown", SPHEROID["WGS84",6378137,298.257223563]], PRIMEM["Greenwich",0], UNIT["degree",0.0174532925199433]], PROJECTION["Transverse_Mercator"], PARAMETER["latitude_of_origin",0], PARAMETER["central_meridian",21], PARAMETER["scale_factor",1], PARAMETER["false_easting",0], PARAMETER["false_northing",0], UNIT["metre",1]]', ...
%     [outImageFileNames{1}(1:end-4) '_GT.tif']));
for i = 1:length(outImageFileNames)
    out = imread(outImageFileNames{i});
%     im = imread(imageFileNames{i});
    gi = geotiffinfo(imageFileNames{i});
    newoutImageFileName = [outImageFileNames{i}(1:end-4) '_GT.tif'];
    %fix no_data
%     out(out>250)=250;
%     mask = repmat(any(out==255, 3), [1, 1, 3]);
%     out(mask) = 2;
%     border = all(im == 0, 3);
%     out(repmat(border, [1, 1, 3])) = 0;
    border = all(out==0, 3);
    out(out==0) = 1;
    out(repmat(border, [1 1 3])) = 0;
    geotiffwrite(newoutImageFileName, out, gi.SpatialRef, ...
        'GeoKeyDirectoryTag', gi.GeoTIFFTags.GeoKeyDirectoryTag);
    fprintf('.');
    %copy and paste into dos win
end
fprintf('\n');


for i = 1:3
    newoutImageFileName = [outImageFileNames{i}(1:end-4) '_GT.tif']; 
    dos(sprintf('gdal_edit.py -a_nodata 0 -a_srs "%s" "%s"\n', ... %
         'PROJCS["unnamed", GEOGCS["unknown", DATUM["unknown", SPHEROID["WGS84",6378137,298.257223563]], PRIMEM["Greenwich",0], UNIT["degree",0.0174532925199433]], PROJECTION["Transverse_Mercator"], PARAMETER["latitude_of_origin",0], PARAMETER["central_meridian",21], PARAMETER["scale_factor",1], PARAMETER["false_easting",0], PARAMETER["false_northing",0], UNIT["metre",1]]', ...
         newoutImageFileName));
end


for i = 4:4
    newoutImageFileName = [outImageFileNames{i}(1:end-4) '_GT.tif'];
    dos(sprintf('gdal_edit.py -a_nodata 0 -a_srs "%s" "%s"\n', ... %-a_nodata 0 
         'PROJCS["unnamed", GEOGCS["unknown", DATUM["unknown", SPHEROID["WGS84",6378137,298.257223563]], PRIMEM["Greenwich",0], UNIT["degree",0.0174532925199433]], PROJECTION["Transverse_Mercator"], PARAMETER["latitude_of_origin",0], PARAMETER["central_meridian",23], PARAMETER["scale_factor",1], PARAMETER["false_easting",0], PARAMETER["false_northing",0], UNIT["metre",1]]', ...
         newoutImageFileName));
end

tmp = imread(newoutImageFileName);
imshow(tmp)

%% -----------------------------------------------------------------------
% Visualise close-up boundary of results

%-------------------------------------------------------------------------
%% Visualise all data by slope

close all; clear all;
dataFileName = 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\DataAllWin5NoBorder2.mat';
prload(dataFileName);

feats = [2 5];
dataAll = changelablist(dataAll, 'Default');
classes = cellstr(getlablist(dataAll));
dataAll = changelablist(dataAll, 'Slope');
slopes = cellstr(getlablist(dataAll));
% slopes(3) = []; %hack out misspelled valleyt

figure;
scatterd2(dataAll(1:10:end, feats), 'legend')
title('All Classes')

dataAll = changelablist(dataAll, 'Default');
figure;
scatterd2(dataAll(1:10:end, feats), 'legend')
title('All Slopes')

figure;
for i = 1:length(classes)
    dataAll = changelablist(dataAll, 'Default');
    classIdx = (getnlab(dataAll) == getclassi2(dataAll, classes{i}));
    subData = dataAll(classIdx,:);
    subData = changelablist(subData, 'Slope');
    h(i) = subplot(2,2,i);
    scatterd2(subData(1:5:end, feats), 'legend')
    title(classes{i})
end
linkaxes(h, 'xy');

err = [];
figure;
for i = 1:length(slopes)
    dataAll = changelablist(dataAll, 'Slope');
    classIdx = (getnlab(dataAll) == getclassi2(dataAll, slopes{i}));
    subData = dataAll(classIdx,:);
    subData = changelablist(subData, 'Default');
    h(i) = subplot(2,3,i);
    scatterd2(subData(1:5:end, feats), 'legend');
    title(slopes{i})
    
    subData = setprior(subData, 0);
    [tr ts] = gendat(subData(:, [2 5]), 0.5);
    if getsize(tr, 3) ==3
        w = tr*qdc;
        c = confmat(ts*w);
        cn = c./repmat(sum(c, 2), 1, size(c, 2));
        err = [err 1-mean(diag(cn))];
    end
end
linkaxes(h, 'xy');
fprintf('Average error for separation by slopes: %f\n', mean(err));



%-------------------------------------------------------------------------
%% Visualise all data by habitat

close all; clear all;
% dataFileNames = {...
%     'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321b_3172_12_0415_rgbn_data.mat';
%     'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321b_3172_12_0419_rgbn_data.mat';
%     'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321D_319_04_0121_rgbn_data.mat';
%     'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3322C_322_02_0056_rgbn_data.mat';
%     };
% dataAll = prdataset();
% 
% for i = 1:length(dataFileNames)
%     load(dataFileNames{i});
%     data = changelablist(data, 'Default');
%     dataAll = [dataAll;data];
% end
prload('G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\DataAllWin5NoBorder2.mat');
dataAll = changelablist(dataAll, 'Default');

feats = [2 5];
classes = cellstr(getlablist(dataAll));
dataAll = changelablist(dataAll, 'Habitat');
habitats = cellstr(getlablist(dataAll));
habitats(3) = []; %hack out misspelled valleyt

figure;
scatterd2(dataAll(1:10:end, feats), 'legend')
title('All Classes')

dataAll = changelablist(dataAll, 'Default');
figure;
scatterd2(dataAll(1:10:end, feats), 'legend')
title('All Habitats')

figure;
for i = 1:length(classes)
    dataAll = changelablist(dataAll, 'Default');
    classIdx = (getnlab(dataAll) == getclassi2(dataAll, classes{i}));
    subData = dataAll(classIdx,:);
    subData = changelablist(subData, 'Habitat');
    h(i) = subplot(2,2,i);
    scatterd2(subData(1:5:end, feats), 'legend')
    title(classes{i})
end
linkaxes(h, 'xy');

err = [];
figure;
for i = 1:length(habitats)
    dataAll = changelablist(dataAll, 'Habitat');
    classIdx = (getnlab(dataAll) == getclassi2(dataAll, habitats{i}));
    subData = dataAll(classIdx,:);
    subData = changelablist(subData, 'Default');
    h(i) = subplot(2,3,i);
    scatterd2(subData(1:5:end, feats), 'legend');
    title(habitats{i})
    
    subData = setprior(subData, 0);
    [tr ts] = gendat(subData(:, [2 5]), 0.5);
    if getsize(tr, 3) ==3
        w = tr*qdc;
        c = confmat(ts*w);
        cn = c./repmat(sum(c, 2), 1, size(c, 2));
        err = [err 1-mean(diag(cn))];
    end
end
linkaxes(h, 'xy');
fprintf('Average error for separation by habitats: %f\n', mean(err));

%-------------------------------------------------------------------------
%% Visualise all data by area 2

clear all;close all;
dataFileNames = {...
    'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321b_3172_12_0415_rgbn_data.mat';
    'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321b_3172_12_0419_rgbn_data.mat';
    'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321D_319_04_0121_rgbn_data.mat';
    'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3322C_322_02_0056_rgbn_data.mat';
    };
dataAll = prdataset();

mainAreas = {'Groenfontein', 'Matjiesvlei', 'Rooiberg', 'GrootKop'};
for i = 1:length(dataFileNames)
    load(dataFileNames{i});
    data = changelablist(data, 'Default');
    data = addlabels(data, mainAreas{i}, 'MainArea');
    dataAll = [dataAll;data];
end

feats = [2 5];
dataAll = changelablist(dataAll, 'Default');
classes = cellstr(getlablist(dataAll));
dataAll = changelablist(dataAll, 'MainArea');
areas = cellstr(getlablist(dataAll));

figure;
scatterd2(dataAll(1:10:end, feats), 'legend')
title('All Classes')

dataAll = changelablist(dataAll, 'Default');
figure;
scatterd2(dataAll(1:10:end, feats), 'legend')
title('All Areas')

figure;
for i = 1:length(classes)
    dataAll = changelablist(dataAll, 'Default');
    classIdx = (getnlab(dataAll) == getclassi2(dataAll, classes{i}));
    subData = dataAll(classIdx,:);
    subData = changelablist(subData, 'MainArea');
    h(i) = subplot(2,2,i);
    scatterd2(subData(1:5:end, feats), 'legend')
    title(classes{i})
end
linkaxes(h, 'xy');

figure;
for i = 1:length(areas)
    dataAll = changelablist(dataAll, 'MainArea');
    classIdx = (getnlab(dataAll) == getclassi2(dataAll, areas{i}));
    subData = dataAll(classIdx,:);
    subData = changelablist(subData, 'Default');
    h(i) = subplot(2,3,i);
    scatterd2(subData(1:5:end, feats), 'legend');
    title(areas{i})

    subData = setprior(subData, 0);
    [tr ts] = gendat(subData(:, feats), 0.5);
    w = tr*qdc;
    c = confmat(ts*w);
    cn = c./repmat(sum(c, 2), 1, size(c, 2));
    err(i) = 1 - mean(diag(cn));
    
%     err(i) = ts*w*testc;
end
linkaxes(h, 'xy');
fprintf('Average error for separation by area: %f\n', mean(err));

%%
load(dataFileNames{1});
data = changelablist(data, 'Default');
classes = cellstr(getlablist(data));
for i = 1:length(classes)+1
    h{i} = [];
    leg{i} = [];
end
cols = {'r','g','b','c','m','y','k'};
icons = {'x','o','<'};

for i = 1:length(dataFileNames)
    load(dataFileNames{i});
    data = changelablist(data, 'Habitat');
    areas_ = cellstr(getlablist(data));
    areas{i} = areas_{1}(1:end-2);
    data = changelablist(data, 'Default');
    dataAll = [dataAll;data];
%     classes = cellstr(getlablist(data));
    for j = 1:length(classes)
        classIdx{j} = (getnlab(data)==getclassi2(data, classes{j}));
        figure(j)
        h_ = plot(+data(classIdx{j}, feats(1)), +data(classIdx{j}, feats(2)), [cols{i} icons{j}]);
        h{j}(end+1) = h_(1);
        leg{j}{end+1} = [areas{i}];
        hold on
        legend(h{j},leg{j});
        title(classes{j})
    end
    figure(length(classes)+1)
    h_ = plot(+data(:, feats(1)), +data(:, feats(2)), [cols{i} icons{1}]);
    h{length(classes)+1}(end+1) = h_(1);
    leg{length(classes)+1}{end+1} = [areas{i}];
    hold on
    
end
figure(length(classes)+1)
legend(h{length(classes)+1},leg{length(classes)+1});

dataAll = changelablist(dataAll, 'Habitat');
figure;
scatterd2(dataAll(1:5:end,feats), 'legend')
title('All areas')

%-------------------------------------------------------------------------
%% Test separation of all data

% prwarning on
warning off
clear all; %close all;
% % % dataFileNames = {...
% % %     'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321b_3172_12_0415_rgbn_data.mat';
% % %     'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321b_3172_12_0419_rgbn_data.mat';
% % %     'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321D_319_04_0121_rgbn_data.mat';
% % %     'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3322C_322_02_0056_rgbn_data.mat';
% % %     };
% % dataAll = prdataset();
% 
% for i = 1:length(dataFileNames)
%     load(dataFileNames{i});
%     dataAll = [dataAll;prdataset(data)];
% end
load('G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\DataAllWin5NoBorder2.mat');
dataAll = changelablist(dataAll, 'Default');
% p = [1 5 1];
% dataAll = setprior(dataAll, p./sum(p));

dataAll = setprior(dataAll, 0);
% data = setprior(data, 0);
% [tr ts] = gendat(data(1:10:end,:), 0.5);
[tr ts] = gendat(dataAll, 0.5);
fl = cellstr(getfeatlab(dataAll));
classes = cellstr(getlablist(dataAll))


[wfs,r] = featselo(tr(1:5:end,:), pca([], .9)*qdc, 4, ts(1:5:end,:));

% scatterdui(tr)
% 
% [wfs,r] = featselo(tr(1:5:end,:), qdc, 2, ts(1:5:end,:));
% getfeatlab(tr*wfs) 
% 
% G           
% rN          
% stdIrRat
% 
% NIR
% irRat
% entropyIrRat
% 
% entropyIrRat qdc
% G           
% rN    
% 
% R           nbc
% NIR         
% NDVI 
% 
% G           nbc
% NDVI  
% 
% G           qdc    
% rN      

feats = [2 4 9];
feats = [13 2 5];
feats = [12 4 1 5];
feats = [5 7 8 2];
feats = [2 5];

w2 = tr(1:5:end, feats)*mogc([], [3 2 2])

w = tr(1:5:end, feats)*(nlfisherm([],2)*qdc)

w = tr(1:50:end, :)*librandomforestc([], 10, 2); %([], 50, 2);

% w2 = tr(:, feats)*qdc
ts(:, feats)*w2*testc
ts*w*testc
% ts(:, feats)*w2*testc
save('G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\All_clfr.mat', 'w', 'feats');

figure;
scatterd2(tr(:, feats), 'legend')
h = plotc(w, 'k-');

c = confmat(ts(:, feats)*w);
cn = c./repmat(sum(c, 2), 1, size(c, 2))
1 - mean(diag(cn))

[CRn] = ReduceConfMat(cn, {[1 3], [2]}, true)
1 - mean(diag(CRn))
%-------------------------------------------------------------------------
%% Test separation of just Spekboom and background

subIdx = getnlab(dataAll) ~= getclassi(dataAll, 'Tree');

subData = dataAll(subIdx, :);

% figure;
scatterdui(subData)
title('All areas')

subData = setprior(subData, 0);
% data = setprior(data, 0);
% [tr ts] = gendat(data(1:10:end,:), 0.5);
[tr ts] = gendat(subData, 0.5);
fl = cellstr(getfeatlab(subData));
classes = cellstr(getlablist(subData))

% scatterdui(tr)
% 
[wfs,r] = featselo(tr(1:5:end,:), qdc, 2, ts(1:5:end,:));
getfeatlab(tr*wfs) 

% bN           qdc
% irRat 
% 
% bN         qdc 
% rN          
% irRat
% 
% bN      nbc      
% rN          
% NDVI  
% 
% irRat    nbc    
% NDVI

feats = [7 5 10];
feats = [7 5 9];

w = tr(:, feats)*qdc
ts(:, feats)*w*testc

c = confmat(ts(:, feats)*w);
cn = c./repmat(sum(c, 2), 1, size(c, 2))
1 - mean(diag(cn))


% [CRn] = ReduceConfMat(cn, {[1 3], [2]}, true)
% 1 - mean(diag(CRn))

%-------------------------------------------------------------------------
%% Test separation of just Spekboom and trees

subIdx = getnlab(dataAll) ~= getclassi(dataAll, 'Background');

subData = dataAll(subIdx, :);

% figure;
scatterdui(subData)
title('All areas')

subData = setprior(subData, 0);
% data = setprior(data, 0);
% [tr ts] = gendat(data(1:10:end,:), 0.5);
[tr ts] = gendat(subData, 0.5);
fl = cellstr(getfeatlab(subData));
classes = cellstr(getlablist(subData))

% scatterdui(tr)
% 
[wfs,r] = featselo(tr(1:2:end,:), qdc, 2, ts(1:2:end,:));
getfeatlab(tr*wfs) 

% G           qdc
% NDVI  
% rN          qdc
% G   

feats = [7 5 10];
feats = [7 5 9];
feats = [2 5];

w = tr(:, feats)*qdc
ts(:, feats)*w*testc

c = confmat(ts(:, feats)*w);
cn = c./repmat(sum(c, 2), 1, size(c, 2))
1 - mean(diag(cn))

figure;
scatterd2(subData(:, feats), 'legend')
plotc(w, 'k-')
% [CRn] = ReduceConfMat(cn, {[1 3], [2]}, true)
% 1 - mean(diag(CRn))


%-------------------------------------------------------------------------
%% OCC approach on all data
% warning off

clear all;close all;
dataFileNames = {...
    'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321b_3172_12_0415_rgbn_data.mat';
    'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321b_3172_12_0419_rgbn_data.mat';
    'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321D_319_04_0121_rgbn_data.mat';
    'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3322C_322_02_0056_rgbn_data.mat';
    };
dataAll = prdataset();

for i = 1:length(dataFileNames)
    load(dataFileNames{i});
    dataAll = [dataAll;data];
end

dataAll = changelablist(dataAll, 'Default');
dataAll = setprior(dataAll, 0);
fl = cellstr(getfeatlab(dataAll));

%Note that the trees and background are not getting divided up with ==
%priors this way
n = min(classsizes(dataAll));
rDataAll = prdataset();
for i = 1:3
    classData = seldat(dataAll, 1);
    classData = gendat(dataAll, n);
    rDataAll = [rDataAll; classData];
end

dataAll = oc_set(rDataAll, 2);
[tr ts] = gendat(dataAll, 0.5);

[wfs,r] = featselo(tr(1:10:end,:), scalem([], 'variance')*gauss_dd, 2, ts(1:10:end,:));
getfeatlab(tr*wfs) 

feats = [2 4 9];
feats = [13 2 5];
feats = [2 5];
% feats = [5 9];
% [tr ts] = gendat(dataAll(:, feats), 0.5);

% train the individual data descriptions and plot them
% the error on the target class:
fracrej = 0.1;

% train the nndd:
% w1 = mcd_gauss_dd(tr, fracrej);
w1 = scalem([], 'variance')*mog_dd([], fracrej, [2 3], 'full');
w1 = tr(:, feats)*w1;
w = w1;
save('G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\All_Occ_clfr.mat', 'w', 'feats');
% w2 = svdd(tr(1:1000:end,:), fracrej, 5);

% and plot the decision boundary:
figure;
scatterd2(tr(:, feats), 'legend')
h = plotc(w1, 'k-');
% h2 = plotc(w2, 'g-');

c = confmat(ts(:, feats)*w1);
cn = c./repmat(sum(c, 2), 1, size(c, 2))
1 - mean(diag(cn))


%--------------------------------------------------------------------------
%% Run classifier on images

close all; clear all;
warning off
fn(1).dataFileName = 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321b_3172_12_0415_rgbn_data.mat';
fn(1).imageFileName = 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321b_3172_12_0415_rgbn_XCALIB.tif';
fn(1).outImageFileName = 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321b_3172_12_0415_rgbn_XCALIB_Out_AllKnnc.tif';
% fn(1).clfrFileName = 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321b_3172_12_0415_clfr.mat';
fn(1).clfrFileName = 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\All_Knnc_clfr.mat';

fn(2).dataFileName = 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321b_3172_12_0419_rgbn_data.mat';
fn(2).imageFileName = 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321b_3172_12_0419_rgbn_XCALIB.tif';
fn(2).outImageFileName = 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321b_3172_12_0419_rgbn_XCALIB_Out_AllKnnc.tif';
% fn(2).clfrFileName = 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321b_3172_12_0419_rgbn_clfr.mat';
fn(2).clfrFileName = 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\All_Knnc_clfr.mat';

fn(3).dataFileName = 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321D_319_04_0121_rgbn_data.mat';
fn(3).imageFileName = 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321D_319_04_0121_rgbn_XCALIB.tif';
fn(3).outImageFileName = 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321D_319_04_0121_rgbn_XCALIB_Out_AllKnnc.tif';
% fn(3).clfrFileName = 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321D_319_04_0121_clfr.mat';
fn(3).clfrFileName = 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\All_Knnc_clfr.mat';

fn(4).dataFileName = 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3322C_322_02_0056_rgbn_data.mat';
fn(4).imageFileName = 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3322C_322_02_0056_rgbn_XCALIB.tif';
fn(4).outImageFileName = 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3322C_322_02_0056_rgbn_XCALIB_Out_AllKnnc.tif';
% fn(4).clfrFileName = 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3322C_322_02_0056_clfr.mat';
fn(4).clfrFileName = 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\All_Knnc_clfr.mat';

% matlabpool
for i = 1:length(fn)
    load(fn(i).clfrFileName);
    blockproc(fn(i).imageFileName, [256 256], @(x)ClassifyIm(x, w, feats), 'Destination', ...
        fn(i).outImageFileName, 'UseParallel', true);    
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
% load(dataFileName);
% data = changelablist(data, 'Default');
% p = [1 2 1];
% data = setprior(data, p./sum(p));
% [tr ts] = gendat(data, 0.5)
% % 
% % feats = [4 9];
% % w = qdc(tr(:, feats))
% 
% load(clfrFileName);
% 
% % confmat(ts(:, feats)*w, 'disagreement')
% 
% figure;
% scatterd2(tr(:, feats), 'legend')
% plotc(w)
% 
% % save('w1', 'w');
% confmat(ts(:, feats)*w, 'disagreement');
% c = confmat(ts(:, feats)*w);
% cn = c./repmat(sum(c, 2), 1, size(c, 2))
% 
% [CRn] = ReduceConfMat(cn, {[1 3], [2]}, true)
% %%
% % matlabpool
% blockproc(imageFileName, [256 256], @(x)ClassifyIm(x, w, feats), 'Destination', ...
%     outImageFileName, 'UseParallel', true);
% % IMWRITE2TIF(IMGDATA,HEADER,IMFILE,DATATYPE)
% 
% % imwrite2tif((out), [], outImageFileName, 'double');
% % save([outImageFileName '.mat'], 'out', '-v7.3')
% % imwrite(out, outImageFileName, 'tiff'); %converts to uint8
% % geotiffinfo(outImageFileName)
% 
% % load(outImageFileName);
% out = imread(outImageFileName);
% 
% [dum outC] = max(+out,[],3);
% outC = outC==2;
% clear dum %out
% 
% % outIm = reshape(+out(:,2), size(mask));
% im = imread(imageFileName);
% %
% figure;
% h1 = subplot(1,3,1);
% imshow(out(:,:,2));
% h2 = subplot(1,3,2);
% imshow(im(:,:,[4 1 2])*16);
% h3 = subplot(1,3,3);
% imshow(im(:,:,[1 2 3])*16);
% linkaxes([h1 h2 h3], 'xy');
% 
