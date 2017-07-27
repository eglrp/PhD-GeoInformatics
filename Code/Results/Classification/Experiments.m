%%
%Visualise prev ground truth
close all; clear all;
cd 'D:\Data\Development\Projects\MSc GeoInformatics\Code\Results\Report1'
rgbFilename = 'F:\MSc GeoInformatics\Data\NGI\Temp\SpekboomThicketCloseupRgb1.png';
cirFilename = 'F:\MSc GeoInformatics\Data\NGI\Temp\SpekboomThicketCloseupCir1.png';
dataFilename = 'F:\MSc GeoInformatics\Data\NGI\Temp\SpekboomThicketCloseupDataset1.mat';

[combinedImRgbIr cirIm] = CombineRgbCir(rgbFilename, cirFilename);

% [data imData] = ExtractFeatures(combinedImRgbIr);

load(dataFilename);

figure;
h1 = subplot(1,2,1);
imshow(cirIm(:,:,[1 2 3]));
hold on
h2 = subplot(1,2,2);
imshow(combinedImRgbIr(:,:,[1 2 3]));
hold on

icons = {'r','g','b'};
for c = 1:max(combinedMask(:))
    maskBound = bwboundaries(combinedMask==c);
    for i = 1:length(maskBound)
        plot(h1,maskBound{i}(:, 2), maskBound{i}(:, 1), icons{c});
        plot(h2,maskBound{i}(:, 2), maskBound{i}(:, 1), icons{c});
    end
% bgMaskBound = bwboundaries(bgMask);
end
linkaxes([h1 h2], 'xy')
%%
%test reading and displaying shapefile ground truth
shapeFileName = 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\Test3.shp';
s = shaperead(shapeFileName);
imageFileName = 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321b_3172_12_0415_rgbn_XCALIB.tif';
gi = geotiffinfo(imageFileName);
im = imread(imageFileName);

[i,j] = map2pix(gi.SpatialRef, s(end).X, s(end).Y);

figure;
imshow(im(:,:,[4 1 2])*16)
hold on;
plot(j, i, 'b', 'LineWidth', 5)

%%
%Test read ground truth and extract features
clear all
shapeFileName = 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321B_3172_12_0415_GroundTruth.shp';
imageFileName = 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321b_3172_12_0415_rgbn_XCALIB.tif';
dataFileName = 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3321b_3172_12_0415_rgbn_data.mat';
s = shaperead(shapeFileName);
gi = geotiffinfo(imageFileName);
im = imread(imageFileName);

[i,j] = map2pix(gi.SpatialRef, [s.X], [s.Y]);
figure;
imshow(im(:,:,[4 1 2])*16)
hold on;
plot(j, i, 'b', 'LineWidth', 5)

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

[data imData] = ExtractFeatures2(im, s, gi.SpatialRef);
save(dataFileName, 'data')

scatterdui(data);

%%
close all; clear all;
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

dataFileName = 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3322C_322_02_0056_rgbn_data.mat';
imageFileName = 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3322C_322_02_0056_rgbn_XCALIB.tif';
outImageFileName = 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3322C_322_02_0056_rgbn_XCALIB_Out.tif';
clfrFileName = 'G:\MSc GeoInformatics\Data\NGI\My Rectified\Ground Truth Images\3322C_322_02_0056_clfr.mat';


load(dataFileName);
data = changelablist(data, 'Default');
p = [1 2 1];
data = setprior(data, p./sum(p));
[tr ts] = gendat(data, 0.5)
% 
% feats = [4 9];
% w = qdc(tr(:, feats))

load(clfrFileName);

% confmat(ts(:, feats)*w, 'disagreement')

figure;
scatterd2(tr(:, feats), 'legend')
plotc(w)

% save('w1', 'w');
confmat(ts(:, feats)*w, 'disagreement');
c = confmat(ts(:, feats)*w);
cn = c./repmat(sum(c, 2), 1, size(c, 2))

[CRn] = ReduceConfMat(cn, {[1 3], [2]}, true)
%%
% matlabpool
blockproc(imageFileName, [256 256], @(x)ClassifyIm(x, w, feats), 'Destination', ...
    outImageFileName, 'UseParallel', true);
% IMWRITE2TIF(IMGDATA,HEADER,IMFILE,DATATYPE)

% imwrite2tif((out), [], outImageFileName, 'double');
% save([outImageFileName '.mat'], 'out', '-v7.3')
% imwrite(out, outImageFileName, 'tiff'); %converts to uint8
% geotiffinfo(outImageFileName)

% load(outImageFileName);
out = imread(outImageFileName);

[dum outC] = max(+out,[],3);
outC = outC==2;
clear dum %out

% outIm = reshape(+out(:,2), size(mask));
im = imread(imageFileName);
%
figure;
h1 = subplot(1,3,1);
imshow(out(:,:,2));
h2 = subplot(1,3,2);
imshow(im(:,:,[4 1 2])*16);
h3 = subplot(1,3,3);
imshow(im(:,:,[1 2 3])*16);
linkaxes([h1 h2 h3], 'xy');

%%
%check classification results against Jan Vlok's ground truth
% shapeFileName = 'D:\Data\Development\Projects\MSc GeoInformatics\Docs\Trimble\GroundTruthCombinedCorrected_TM21_WGS84.shp';
shapeFileName = 'D:\Data\Development\Projects\MSc GeoInformatics\Docs\Trimble\GroundTruthCombinedCorrected_TM23_WGS84.shp';
s = shaperead(shapeFileName);
gi = geotiffinfo(imageFileName);
% im = imread(imageFileName);
% out = imread(outImageFileName);

for i = 1:length(s)
    [row,col] = map2pix(gi.SpatialRef, [s(i).X], [s(i).Y]);
    if (min(row)>0 && max(row)<=size(im,1) && min(col)>0 && max(col)<=size(im,2))
        rowIdx = floor(min(row)):ceil(max(row));
        colIdx = floor(min(col)):ceil(max(col));

        idx = find(isnan(s(i).X));

%         idx = 1:idx(1)-1;
        mask = false(size(im,1), size(im,2));
        startIdx = 1;
        for j = 1:length(idx)
            polyIdx = startIdx:idx(j)-1;
            mask = mask | poly2mask(col(polyIdx), row(polyIdx), size(im,1), size(im,2));
            startIdx = idx(j)+1;
        end
        fprintf('%s %d-%d, %f\n', s(i).Name, s(i).CoverMin, s(i).CoverMax, 100*sum(outC(mask))/sum(mask(:)));

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

%%
