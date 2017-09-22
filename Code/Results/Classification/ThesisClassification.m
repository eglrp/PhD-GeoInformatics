close all hidden; clear all;
load('D:\Data\Development\Projects\PhD GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\DataAllWin5NoBorder2.mat')
dataAll = changelablist(dataAll, 'Default')
dataAll = setprior(dataAll, 0);
fl = cellstr(getfeatlab(dataAll))

global feats
feats = [9 15 20 23 7 6]; %ranked cluster
% feats = [16 23 21 25 38 20]; %wfso = featself(subData, libsvc([], proxm([], 'r', 2), 20), 6);
% feats = +wfso
fl(feats)
% feats = [9 15 23 7 6];
% subData = gendat(dataAll, [2000 2000 2000]);
% subData = changelablist(subData, 'Default');
% subData = setprior(subData, 0);

cs = classsizes(dataAll);
cs(1) = cs(2);
subData = gendat(dataAll, cs);
subData = changelablist(subData, 'Default');
subData = setprior(subData, 0);

%% ----------------------------- Decision tree
% d = gendat(dataAll, [2000 2000 2000]);
% d = setprior(d, 0);
d = subData;
tic
p = [1 1 1];
[err, cerr, nlabOut] = prcrossval(d(:, feats), opencvdtreec([], 12, {'Priors', p./sum(p), 'MaxDepth', 12, 'Use1seRule', false, ...
        'UseSurrogates', false, 'CVFolds', 5, 'TruncatePrunedTree', true, 'MinSampleCount', min(classsizes(subData))/100}), 10);
toc
[cc ccr] = ClfrPerfMeas(d, nlabOut);

%% ----------------------------- KNN
% [d, dts] = gendat(dataAll, [2000 2000 2000]);
d = subData;

[err, cerr, nlabOut, stds, r] = prcrossval(d(:, feats), scalem([], 'variance')*knnc([], 5), 10);

ClfrPerfMeas(d, nlabOut);

%% ----------------------------- Random forest on selected features
d = subData;

[err, cerr, nlabOut, stds, r] = prcrossval(d(:, feats), opencvrtreec([], [], {'Priors', [1 2 1]/4, ...
     'MaxNumOfTreesInTheForest', 5, 'NActiveVars', 4, 'CalcVarImportance', false, 'MaxDepth', 10, 'ForestAccuracy', 0.025}), 10);

ClfrPerfMeas(d, nlabOut);

%% ----------------------------- Random forest on full feature set
d = subData;
p = [1 1.5 1];
[err, cerr, nlabOut, stds, r] = prcrossval(d, opencvrtreec([], [], {'Priors', p./sum(p), ...
     'MaxNumOfTreesInTheForest', 50, 'NActiveVars', 4, 'CalcVarImportance', false, 'MaxDepth', 10, 'ForestAccuracy', 0.025}), 2);

ClfrPerfMeas(d, nlabOut);

%% ----------------------------- libsvc on selected features
% d = gendat(dataAll, [2000 2000 2000]);
% d = setprior(d, [0.25 0.5 0.25]);
d = subData;

[err, cerr, nlabOut, stds, r] = prcrossval(d(:, feats), scalem([], 'variance')*libsvc([], proxm([], 'r', 2), 20), 10);

ClfrPerfMeas(d, nlabOut);


%% ----------------------------- opencvsvc on selected features
d = subData;

[err, cerr, nlabOut, stds, r] = prcrossval(d(:, feats), scalem([], 'variance')*opencvsvc([], [], ...
    {'SVMType', 'C_SVC', 'KernelType', 'RBF', 'Gamma', 25, 'C', 1, 'ClassWeights', double([1; 1; 1])}), 10);

ClfrPerfMeas(d, nlabOut);
%% ----------------------------- qdc on selected features
d = subData;
p = [1 1 1];
d = setprior(d, p./sum(p))
tic
[err, cerr, nlabOut, stds, r] = prcrossval(d(:, feats), qdc, 10);
toc

ClfrPerfMeas(d, nlabOut);

%% ----------------------------- test externally gen output tifs against field GT
%% Validation - get results for all jan vlok gt

close all; 
imageFileNames = {...
    'D:\Data\Development\Projects\PhD GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321b_3172_12_0415_rgbn_XCALIB.tif';...
    'D:\Data\Development\Projects\PhD GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321b_3172_12_0419_rgbn_XCALIB.tif';...
    'D:\Data\Development\Projects\PhD GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321D_319_04_0121_rgbn_XCALIB.tif';...
    'D:\Data\Development\Projects\PhD GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3322C_322_02_0056_rgbn_XCALIB.tif';...
    };
outImageFileNames = {...
    'F:\MSc GeoInformatics\Data\NGI\Classification\3321b_3172_12_0415_rgbn_XCALIB_Out_DTree3.tif';...
    'F:\MSc GeoInformatics\Data\NGI\Classification\3321b_3172_12_0419_rgbn_XCALIB_Out_DTree3.tif';...
    'F:\MSc GeoInformatics\Data\NGI\Classification\3321D_319_04_0121_rgbn_XCALIB_Out_DTree3.tif';...
    'F:\MSc GeoInformatics\Data\NGI\Classification\3322C_322_02_0056_rgbn_XCALIB_Out_DTree3.tif';...
    };

outImageFileNames = {...
    'F:\MSc GeoInformatics\Data\NGI\Classification\3321b_3172_12_0415_rgbn_XCALIB_Out_KNearest2.tif';...
    'F:\MSc GeoInformatics\Data\NGI\Classification\3321b_3172_12_0419_rgbn_XCALIB_Out_KNearest2.tif';...
    'F:\MSc GeoInformatics\Data\NGI\Classification\3321D_319_04_0121_rgbn_XCALIB_Out_KNearest2.tif';...
    'F:\MSc GeoInformatics\Data\NGI\Classification\3322C_322_02_0056_rgbn_XCALIB_Out_KNearest2.tif';...
    };

outImageFileNames = {...
    'F:\MSc GeoInformatics\Data\NGI\Classification\3321b_3172_12_0415_rgbn_XCALIB_Out_SvmOpenCv2.tif';...
    'F:\MSc GeoInformatics\Data\NGI\Classification\3321b_3172_12_0419_rgbn_XCALIB_Out_SvmOpenCv2.tif';...
    'F:\MSc GeoInformatics\Data\NGI\Classification\3321D_319_04_0121_rgbn_XCALIB_Out_SvmOpenCv2.tif';...
    'F:\MSc GeoInformatics\Data\NGI\Classification\3322C_322_02_0056_rgbn_XCALIB_Out_SvmOpenCv2.tif';...
    };

outImageFileNames = {...
    'F:\MSc GeoInformatics\Data\NGI\Classification\3321b_3172_12_0415_rgbn_XCALIB_Out_RTreesSmall2.tif';...
    'F:\MSc GeoInformatics\Data\NGI\Classification\3321b_3172_12_0419_rgbn_XCALIB_Out_RTreesSmall2.tif';...
    'F:\MSc GeoInformatics\Data\NGI\Classification\3321D_319_04_0121_rgbn_XCALIB_Out_RTreesSmall2.tif';...
    'F:\MSc GeoInformatics\Data\NGI\Classification\3322C_322_02_0056_rgbn_XCALIB_Out_RTreesSmall2.tif';...
    };

outImageFileNames = {...
    'F:\MSc GeoInformatics\Data\NGI\Classification\3321b_3172_12_0415_rgbn_XCALIB_Out_NormalBayes2.tif';...
    'F:\MSc GeoInformatics\Data\NGI\Classification\3321b_3172_12_0419_rgbn_XCALIB_Out_NormalBayes2.tif';...
    'F:\MSc GeoInformatics\Data\NGI\Classification\3321D_319_04_0121_rgbn_XCALIB_Out_NormalBayes2.tif';...
    'F:\MSc GeoInformatics\Data\NGI\Classification\3322C_322_02_0056_rgbn_XCALIB_Out_NormalBayes2.tif';...
    };

gtFileNames = {...
    'C:\Data\Development\Projects\MSc GeoInformatics\Docs\Trimble\GroundTruthCombinedCorrected_TM21_WGS84.shp';...
    'C:\Data\Development\Projects\MSc GeoInformatics\Docs\Trimble\GroundTruthCombinedCorrected_TM21_WGS84.shp';...
    'C:\Data\Development\Projects\MSc GeoInformatics\Docs\Trimble\GroundTruthCombinedCorrected_TM21_WGS84.shp';...
    'C:\Data\Development\Projects\MSc GeoInformatics\Docs\Trimble\GroundTruthCombinedCorrected_TM23_WGS84.shp';...
    };
% 
% clfrFileNames = {...
%     'D:\Data\Development\Projects\PhD GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\All_RandomForest_clfr.mat';...
%     'D:\Data\Development\Projects\PhD GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\All_RandomForest_clfr.mat';...
%     'D:\Data\Development\Projects\PhD GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\All_RandomForest_clfr.mat';...
%     'D:\Data\Development\Projects\PhD GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\All_RandomForest_clfr.mat';...
%     };
%%
res = [];
count = 1;
for g = 1:length(gtFileNames)
    s = shaperead(gtFileNames{g});
    gi = geotiffinfo(imageFileNames{g});
%     im = imread(imageFileNames{g});
    out = imread(outImageFileNames{g});
%     [dum outC] = max(+out, [], 3);
%     outC = out==2; %1;NB
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
            
            tl = {'Severe', 'Moderate', 'Pristine'};
            tn = find(strcmpi(s(i).Transforma, tl));
            
            res(count).Transform = tn;
%             res(count).Transform = s(i).Transforma;
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

figure;
plot(abs([res.gtCover]-[res.clfCover]), [res.Transform], 'x');
xlabel('err')
ylabel('xform')

%% ----------------------------- test matlab clfr against field GT
% extracts only features in 'feats' only on gt polygons
% runs either clfr in clfrFileNames or 'w' if clfrFileNames empty
% probably wont give same res as when run in OpenCV
% also dont forget postproc morphology done in OpenCV

% TO DO: check if we get matlab dtree res matching with ext gen opencv results
% then we can do all comparisons internal to matlab
% Otherwise check if we export a full data set and tr with this in opencv
% if it improves on the [2000 2000 2000] set results.  If not just use the
% exisint *_Out_*.tif files for thesis results.  If major improvement then
% start again...

clear all;
close all hidden; 

imageFileNames = {...
    'D:\Data\Development\Projects\PhD GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321b_3172_12_0415_rgbn_XCALIB.tif';...
    'D:\Data\Development\Projects\PhD GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321b_3172_12_0419_rgbn_XCALIB.tif';...
    'D:\Data\Development\Projects\PhD GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321D_319_04_0121_rgbn_XCALIB.tif';...
    'D:\Data\Development\Projects\PhD GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3322C_322_02_0056_rgbn_XCALIB.tif';...
    };

gtFileNames = {...
    'C:\Data\Development\Projects\PhD GeoInformatics\Docs\Trimble\GroundTruthCombinedCorrected_TM21_WGS84.shp';...
    'C:\Data\Development\Projects\PhD GeoInformatics\Docs\Trimble\GroundTruthCombinedCorrected_TM21_WGS84.shp';...
    'C:\Data\Development\Projects\PhD GeoInformatics\Docs\Trimble\GroundTruthCombinedCorrected_TM21_WGS84.shp';...
    'C:\Data\Development\Projects\PhD GeoInformatics\Docs\Trimble\GroundTruthCombinedCorrected_TM23_WGS84.shp';...
    };

% clfrFileNames = {...
%     'D:\Data\Development\Projects\PhD GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\All_Knnc_clfr.mat';...
%     'D:\Data\Development\Projects\PhD GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\All_Knnc_clfr.mat';...
%     'D:\Data\Development\Projects\PhD GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\All_Knnc_clfr.mat';...
%     'D:\Data\Development\Projects\PhD GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\All_Knnc_clfr.mat';...
%     };
global w feats
clfrFileNames = {... %empty to use workspace vars i.e. feats and w
    '';...
    '';...
    '';...
    '';...
    };

%%
load('D:\Data\Development\Projects\PhD GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\DataAllWin5NoBorder2.mat')
dataAll = changelablist(dataAll, 'Default')
dataAll = setprior(dataAll, 0);
fl = cellstr(getfeatlab(dataAll))

global feats
feats = [9 15 20 23 7 6]; %ranked cluster
% feats = [16 23 21 25 38 20]; %wfso = featself(subData, libsvc([], proxm([], 'r', 2), 20), 6);
% feats = +wfso
fl(feats)
% feats = [9 15 23 7 6];
if false
    subData = gendat(dataAll, [2000 2000 2000]);
    subData = changelablist(subData, 'Default');
    subData = setprior(subData, 0);
else
    cs = classsizes(dataAll);
    cs(1) = cs(2);
    subData = gendat(dataAll, cs);
    subData = changelablist(subData, 'Default');
    subData = setprior(subData, 0);
end

%%
d = subData;

p = [1 1 1];
% [err, cerr, nlabOut] = prcrossval(d(:, feats), opencvdtreec([], 12, {'Priors', p./sum(p), 'MaxDepth', 12, 'Use1seRule', false, ...
%         'UseSurrogates', false, 'CVFolds', 5, 'TruncatePrunedTree', true, 'MinSampleCount', min(classsizes(subData))/100}), 10);
% 
% [err, cerr, nlabOut, stds, r] = prcrossval(d(:, feats), scalem([], 'variance')*knnc([], 5), 2);
% 
% [err, cerr, nlabOut, stds, r] = prcrossval(d(:, feats), opencvrtreec([], [], {'Priors', [1 2 1]/4, ...
%      'MaxNumOfTreesInTheForest', 5, 'NActiveVars', 4, 'CalcVarImportance', false, 'MaxDepth', 10, 'ForestAccuracy', 0.025}), 10);
% 
% [err, cerr, nlabOut, stds, r] = prcrossval(d(:, feats), scalem([], 'variance')*libsvc([], proxm([], 'r', 2), 20), 10);
% 
% [err, cerr, nlabOut, stds, r] = prcrossval(d(:, feats), scalem([], 'variance')*knnc([], 5), 2);
% 
% 
% [err, cerr, nlabOut, stds, r] = prcrossval(d(:, feats), scalem([], 'variance')*opencvsvc([], [], ...
%     {'SVMType', 'C_SVC', 'KernelType', 'RBF', 'Gamma', 10, 'C', 1, 'ClassWeights', double([2; 1; 5])}), 2);
p = [1 1 1];
feats = 1:size(subData)
w = d(:, feats)*opencvdtreec([], 12, {'Priors', p./sum(p), 'MaxDepth', 12, 'Use1seRule', false, ...
        'UseSurrogates', false, 'CVFolds', 5, 'TruncatePrunedTree', true, 'MinSampleCount', min(classsizes(subData))/100})

w = subData(:, feats)*(scalem([], 'variance')*knnc([], 5));
w = subData(:, feats)*opencvrtreec([], [], {'Priors', [1 2 1]/4, ...
      'MaxNumOfTreesInTheForest', 5, 'NActiveVars', 4, 'CalcVarImportance', false, 'MaxDepth', 10, 'ForestAccuracy', 0.025})

p = [1 1 1];
w = subData(:, feats)*(scalem([], 'variance')*opencvsvc([], [],...
    {'SVMType', 'C_SVC', 'KernelType', 'RBF', 'Gamma', 25, 'C', .1, 'ClassWeights', double(p)./sum(p)}))

% w = subData(:, feats)*(scalem([], 'variance')*libsvc([], proxm([], 'r', 2), 20));

w = subData(:, feats)*qdc

c = confmat(subData(:, feats)*w);
cn = c./repmat(sum(c,2),1,3)

  %%
% ClfrPerfMeas(d, nlabOut);
close all hidden;
resC = {};
res = [];
for g = 1:length(gtFileNames)
    resC{g} = ValidateClfrAgainstJvGroundTruth(imageFileNames{g}, gtFileNames{g}, clfrFileNames{g});
end

res = [resC{:}];

for i = 1:length(res)
      fprintf('%s %f %f\n', res(i).name, res(i).gtCover, res(i).clfCover);
end

fprintf('Canopy cover error Mean(abs) Std(abs): %f (%f)\n', mean(abs([res.gtCover]-[res.clfCover])), std(abs([res.gtCover]-[res.clfCover])));
fprintf('Canopy cover error Median Mad: %f (%f)\n', median([res.gtCover]-[res.clfCover]), mad([res.gtCover]-[res.clfCover]));

%% --------------------------------- Run OpenCV exe on ground truth tiffs

close all; 
imageFileNames = {...
    'D:\Data\Development\Projects\PhD GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321b_3172_12_0415_rgbn_XCALIB.tif';...
    'D:\Data\Development\Projects\PhD GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321b_3172_12_0419_rgbn_XCALIB.tif';...
    'D:\Data\Development\Projects\PhD GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321D_319_04_0121_rgbn_XCALIB.tif';...
    'D:\Data\Development\Projects\PhD GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3322C_322_02_0056_rgbn_XCALIB.tif';...
    };

%%
    cd('D:\Data\Development\Projects\MSc GeoInformatics\Code\SpekboomClassifier\x64\Release');
    
    outImageFileNames = {...
        'F:\MSc GeoInformatics\Data\NGI\Classification\3321b_3172_12_0415_rgbn_XCALIB_Out_KNearest2.tif';...
        'F:\MSc GeoInformatics\Data\NGI\Classification\3321b_3172_12_0419_rgbn_XCALIB_Out_KNearest2.tif';...
        'F:\MSc GeoInformatics\Data\NGI\Classification\3321D_319_04_0121_rgbn_XCALIB_Out_KNearest2.tif';...
        'F:\MSc GeoInformatics\Data\NGI\Classification\3322C_322_02_0056_rgbn_XCALIB_Out_KNearest2.tif';...
        };
    
    for i = 1:length(imageFileNames)
        fprintf('SpekboomClassifier.exe "%s" "%s" -o -t 1\n', imageFileNames{i}, outImageFileNames{i});
    end
    
    outImageFileNames = {...
        'F:\MSc GeoInformatics\Data\NGI\Classification\3321b_3172_12_0415_rgbn_XCALIB_Out_SvmOpenCv2.tif';...
        'F:\MSc GeoInformatics\Data\NGI\Classification\3321b_3172_12_0419_rgbn_XCALIB_Out_SvmOpenCv2.tif';...
        'F:\MSc GeoInformatics\Data\NGI\Classification\3321D_319_04_0121_rgbn_XCALIB_Out_SvmOpenCv2.tif';...
        'F:\MSc GeoInformatics\Data\NGI\Classification\3322C_322_02_0056_rgbn_XCALIB_Out_SvmOpenCv2.tif';...
        };
    
    for i = 1:length(imageFileNames)
        fprintf('SpekboomClassifier.exe "%s" "%s" -o -t 2\n', imageFileNames{i}, outImageFileNames{i});
    end

    outImageFileNames = {...
        'F:\MSc GeoInformatics\Data\NGI\Classification\3321b_3172_12_0415_rgbn_XCALIB_Out_RTreesSmall2.tif';...
        'F:\MSc GeoInformatics\Data\NGI\Classification\3321b_3172_12_0419_rgbn_XCALIB_Out_RTreesSmall2.tif';...
        'F:\MSc GeoInformatics\Data\NGI\Classification\3321D_319_04_0121_rgbn_XCALIB_Out_RTreesSmall2.tif';...
        'F:\MSc GeoInformatics\Data\NGI\Classification\3322C_322_02_0056_rgbn_XCALIB_Out_RTreesSmall2.tif';...
        };
    
    for i = 1:length(imageFileNames)
        fprintf('SpekboomClassifier.exe "%s" "%s" -o -t 0\n', imageFileNames{i}, outImageFileNames{i});
    end
    
    outImageFileNames = {...
        'F:\MSc GeoInformatics\Data\NGI\Classification\3321b_3172_12_0415_rgbn_XCALIB_Out_NormalBayes2.tif';...
        'F:\MSc GeoInformatics\Data\NGI\Classification\3321b_3172_12_0419_rgbn_XCALIB_Out_NormalBayes2.tif';...
        'F:\MSc GeoInformatics\Data\NGI\Classification\3321D_319_04_0121_rgbn_XCALIB_Out_NormalBayes2.tif';...
        'F:\MSc GeoInformatics\Data\NGI\Classification\3322C_322_02_0056_rgbn_XCALIB_Out_NormalBayes2.tif';...
        };
    
    for i = 1:length(imageFileNames)
        fprintf('SpekboomClassifier.exe "%s" "%s" -o -t 3\n', imageFileNames{i}, outImageFileNames{i});
    end
    
%% -------------------------------------- Compare window sizes
close hidden; clear all;
dataFileNames = {...
    'D:\Data\Development\Projects\PhD GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\DataAllWin3.mat';...
    'D:\Data\Development\Projects\PhD GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\DataAllWin5.mat';...
    'D:\Data\Development\Projects\PhD GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\DataAllWin7.mat';...
    'D:\Data\Development\Projects\PhD GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\DataAllWin9.mat';...
    'D:\Data\Development\Projects\PhD GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\DataAllWin11.mat';...
    'D:\Data\Development\Projects\PhD GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\DataAllWin13.mat';...
    };

featIdx = 15;
feats = [6 9 featIdx];
p = [1 1 1];
for i = 1:length(dataFileNames)
    prload(dataFileNames{i});
    dataAll = changelablist(dataAll, 'Default');
    dataAll = setprior(dataAll, 0);
    cs = classsizes(dataAll);
    cs(1) = cs(2);
    subData = gendat(dataAll, cs);
 
    w = opencvdtreec([], 12, {'Priors', p./sum(p), 'MaxDepth', 12, 'Use1seRule', false, ...
        'UseSurrogates', false, 'CVFolds', 5, 'TruncatePrunedTree', true, 'MinSampleCount', min(classsizes(subData))/100});
%     w = naivebc;
    [err(i), cerr_] = prcrossval(subData(:, feats), w, 10);
% , 10, 1, testc);
%     cerr(i) = mean(cerr_)
    [tr, ts] = gendat(subData, 0.5);
    w = tr(:, feats)*naivebc;
    err2(i) = ts(:, feats)*w*testc;
    fl = cellstr(getfeatlab(dataAll));
    fl(featIdx)
end

figure
plot([3 5 7 9 11 13], err)
hold all
plot([3 5 7 9 11 13], err2)
legend({'err', 'err2'});
t = char(fl(feats))';
title(t(:)')
%Note
% - NB Results here are v data sensitive i.e. repetitive runs do not yield
% same results

%% --------------------------------- See how stratified (on area) sampling changes rooiberg and grootkop performance
close all; clear all;
prload('D:\Data\Development\Projects\PhD GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\DataAllWin5NoBorder2.mat');
dataAll = changelablist(dataAll, 'Default');
dataAll = setprior(dataAll, 0);
feats = [9 15 20 23 7 6]; %ranked cluster

nl = getnlab(dataAll);
ll = cellstr(getlablist(dataAll));

%% Compare separation by custom match list on XX
% compareList = {'Valley', 'Arid'};
% compareBy = 'Habitat';
compareList = {'GroenFontein', 'MatjiesVlei', 'RooiBerg', 'GrootKop'};
compareBy = 'Area';

dataAll = changelablist(dataAll, compareBy);
% sepFields = cellstr(getlablist(dataAll));
sepLabels = cellstr(getlabels(dataAll));
% sepFields(3) = []; %hack out misspelled valleyt

classLabels = {};

for i = 1:length(compareList)
    dataAll = changelablist(dataAll, compareBy);
%     classIdx = (getnlab(dataAll) == getclassi2(dataAll, classes{i}));
    classIdx = strmatch(lower(compareList{i}), lower(sepLabels));
    classLabels(classIdx) = compareList(i);
end

dataAll = addlabels(dataAll, char(classLabels'), 'Area2');
dataAll = changelablist(dataAll, 'Area2')

nlArea = getnlab(dataAll);
llArea = cellstr(getlablist(dataAll));
for i = 1:length(llArea)
    for j = 1:length(ll)
        fprintf('%s - %s: %d %f\n', llArea{i}, ll{j}, sum(nlArea==i & nl==j), sum(nlArea==i & nl==j)/size(dataAll,1));
    end
end

%NOTES
% - Actually the data is quite well stratified already so the reason for the
% poor performance in GrootKop and Rooiberg must be something else

%%

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
% clfr.save('D:\Data\Development\Projects\PhD GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\configDTree.yaml');

clfrSvm = +w;
clfrSvm = +clfrSvm{2};
% clfrSvm.save('D:\Data\Development\Projects\PhD GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\configSvm.yaml');

% save('D:\Data\Development\Projects\PhD GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\All_Svc_clfr3.mat', 'w', 'feats')
% save('D:\Data\Development\Projects\PhD GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\All_RandomForest_clfr_3class.mat', 'w', 'feats')
% save('D:\Data\Development\Projects\PhD GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\All_OpenCvRTrees_clfr_3class.mat', 'w', 'feats')
% save('D:\Data\Development\Projects\PhD GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\All_OpenCvDTree_clfr_3class.mat', 'w', 'feats')
% save('D:\Data\Development\Projects\PhD GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\All_OpenCvSvc_clfr_3class.mat', 'w', 'feats')
% save('D:\Data\Development\Projects\PhD GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\All_OpenCvDTree_clfr_3class.mat', 'w', 'feats', 'subData')

s = [];
s.tr = single(+subData(:, feats));
s.trLab = single(+getnlab(subData));

%copy prev file contents before exec this line
% cv.FileStorage('D:\Data\Development\Projects\PhD GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\trDataDTree.yaml', s);
%%
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

% clfr.save('D:\Data\Development\Projects\PhD GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\configDtree.yaml');

