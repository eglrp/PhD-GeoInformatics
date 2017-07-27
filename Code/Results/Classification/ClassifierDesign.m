%% Feature Extraction
clear all; close all;
% shapeFileName = 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321B_3172_12_0415_GroundTruth.shp';
% imageFileName = 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321b_3172_12_0415_rgbn_XCALIB.tif';
% dataFileName = 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321b_3172_12_0415_rgbn_data.mat';
% clfrFileName = 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321b_3172_12_0415_clfr.mat';

% shapeFileName = 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321B_3172_12_0419_RGBN_GroundTruth.shp';
% imageFileName = 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321b_3172_12_0419_rgbn_XCALIB.tif';
% dataFileName = 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321b_3172_12_0419_rgbn_data.mat';
% clfrFileName = 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321b_3172_12_0419_clfr.mat';

shapeFileName = 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321B_319_04_0121_GroundTruth.shp';
imageFileName = 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321D_319_04_0121_rgbn_XCALIB.tif';
dataFileName = 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321D_319_04_0121_rgbn_data.mat';
clfrFileName = 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321D_319_04_0121_clfr.mat';

% shapeFileName = 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3322C_322_02_0056_GroundTruth.shp';
% imageFileName = 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3322C_322_02_0056_rgbn_XCALIB.tif';
% dataFileName = 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3322C_322_02_0056_rgbn_data.mat';
% clfrFileName = 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3322C_322_02_0056_clfr.mat';


%%
s = shaperead(shapeFileName);
gi = geotiffinfo(imageFileName);
im = imread(imageFileName);

if 1
[i,j] = map2pix(gi.SpatialRef, [s.X], [s.Y]);
% figure;
% imshow(im(:,:,[4 1 2])*16)
% hold on;
% plot(j, i, 'b', 'LineWidth', 5)

combinedMask = zeros(size(im,1), size(im,2));
% bgMask = false(size(im,1), size(im,2));
for i = 1:length(s)
    idx = find(isnan(s(i).X));
    idx = 1:idx(1)-1;
    if (strcmpi(s(i).Class, 'Spekboom'))
        combinedMask = combinedMask | 1*poly2mask(s(i).X(idx), s(i).Y(idx), size(im,1), size(im,2));
    elseif (strcmpi(s(i).Class, 'Background'))
        combinedMask = combinedMask | 2*poly2mask(s(i).X(idx), s(i).Y(idx), size(im,1), size(im,2));
    else
        combinedMask = combinedMask | 3*poly2mask(s(i).X(idx), s(i).Y(idx), size(im,1), size(im,2));
    end
    fprintf('.');
end
fprintf('\n');
end

[data imData] = ExtractFeatures2(im, s, gi.SpatialRef);
save(dataFileName, 'data')

data = changelablist(data, 'Default')
scatterdui(data);

%%
% Visualise features
load(dataFileName)
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
            h_ = plot(+data(idx, feats(1)), +data(idx, feats(2)), [cols{i} icons{j}]);
            h(end+1) = h_(1)
            leg{end+1} = [areas{i} '-' classes{j}];
            hold on
        end
    end
end
legend(h,leg);


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


%%
% Feature selection

load(dataFileName);
data = changelablist(data, 'Default')
p = [1 50 1];
data = setprior(data, p./sum(p));
% data = setprior(data, 0);
% [tr ts] = gendat(data(1:10:end,:), 0.5);
[tr ts] = gendat(data, 0.5);
fl = cellstr(getfeatlab(data));

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


%% -----------------------------------------------------------------------
%Validation - get results for all jan vlok gt
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
    'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321b_3172_12_0415_rgbn_XCALIB_Out_AllOcc.tif';...
    'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321b_3172_12_0419_rgbn_XCALIB_Out_AllOcc.tif';...
    'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321D_319_04_0121_rgbn_XCALIB_Out_AllOcc.tif';...
    'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3322C_322_02_0056_rgbn_XCALIB_Out_AllOcc.tif';...
    };

gtFileNames = {...
    'D:\Data\Development\Projects\MSc GeoInformatics\Docs\Trimble\GroundTruthCombinedCorrected_TM21_WGS84.shp';...
    'D:\Data\Development\Projects\MSc GeoInformatics\Docs\Trimble\GroundTruthCombinedCorrected_TM21_WGS84.shp';...
    'D:\Data\Development\Projects\MSc GeoInformatics\Docs\Trimble\GroundTruthCombinedCorrected_TM21_WGS84.shp';...
    'D:\Data\Development\Projects\MSc GeoInformatics\Docs\Trimble\GroundTruthCombinedCorrected_TM23_WGS84.shp';...
    };

%%
% shapeFileName = 'D:\Data\Development\Projects\MSc GeoInformatics\Docs\Trimble\GroundTruthCombinedCorrected_TM21_WGS84.shp';

count = 1;
for g = 1:length(gtFileNames)
    s = shaperead(gtFileNames{g});
    gi = geotiffinfo(imageFileNames{g});
%     im = imread(imageFileNames{g});
    out = imread(outImageFileNames{g});
    [dum outC] = max(+out, [], 3);
    outC = outC==1; %2;NB
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
%% Visualise all data by habitat

close all; clear all;
dataFileNames = {...
    'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321b_3172_12_0415_rgbn_data.mat';
    'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321b_3172_12_0419_rgbn_data.mat';
    'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321D_319_04_0121_rgbn_data.mat';
    'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3322C_322_02_0056_rgbn_data.mat';
    };
dataAll = dataset();

for i = 1:length(dataFileNames)
    load(dataFileNames{i});
    data = changelablist(data, 'Default');
    dataAll = [dataAll;data];
end

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
dataAll = dataset();

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
dataFileNames = {...
    'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321b_3172_12_0415_rgbn_data.mat';
    'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321b_3172_12_0419_rgbn_data.mat';
    'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321D_319_04_0121_rgbn_data.mat';
    'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3322C_322_02_0056_rgbn_data.mat';
    };
dataAll = prdataset();

for i = 1:length(dataFileNames)
    load(dataFileNames{i});
    dataAll = [dataAll;prdataset(data)];
end
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

w = tr(1:50:end, :)*randomforestc([], 10, 2); %([], 50, 2);

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
dataAll = dataset();

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
rDataAll = dataset();
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
fn(1).outImageFileName = 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321b_3172_12_0415_rgbn_XCALIB_Out_AllQdc.tif';
% fn(1).clfrFileName = 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321b_3172_12_0415_clfr.mat';
fn(1).clfrFileName = 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\All_clfr.mat';

fn(2).dataFileName = 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321b_3172_12_0419_rgbn_data.mat';
fn(2).imageFileName = 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321b_3172_12_0419_rgbn_XCALIB.tif';
fn(2).outImageFileName = 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321b_3172_12_0419_rgbn_XCALIB_Out_AllQdc.tif';
% fn(2).clfrFileName = 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321b_3172_12_0419_rgbn_clfr.mat';
fn(2).clfrFileName = 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\All_clfr.mat';

fn(3).dataFileName = 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321D_319_04_0121_rgbn_data.mat';
fn(3).imageFileName = 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321D_319_04_0121_rgbn_XCALIB.tif';
fn(3).outImageFileName = 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321D_319_04_0121_rgbn_XCALIB_Out_AllQdc.tif';
% fn(3).clfrFileName = 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321D_319_04_0121_clfr.mat';
fn(3).clfrFileName = 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\All_clfr.mat';

fn(4).dataFileName = 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3322C_322_02_0056_rgbn_data.mat';
fn(4).imageFileName = 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3322C_322_02_0056_rgbn_XCALIB.tif';
fn(4).outImageFileName = 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3322C_322_02_0056_rgbn_XCALIB_Out_AllQdc.tif';
% fn(4).clfrFileName = 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3322C_322_02_0056_clfr.mat';
fn(4).clfrFileName = 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\All_clfr.mat';

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
