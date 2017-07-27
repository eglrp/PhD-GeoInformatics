function res = CrossCalib(refFileName, srcFileName, varargin)

refR = [];
calibR = [];
winSize = 1;
res = [];

ModifyDefaultArgs(varargin);

if ischar(refFileName)
    [refIm refR] = geotiffread(refFileName);
else
    refIm = refFileName;
end

if ischar(srcFileName)
    %downsample calibIm with GDAL
    srcInfo = geotiffinfo(srcFileName);
    refPath = fileparts(refFileName);
    [srcPath srcName] = fileparts(srcFileName);
    srcDsFileName = sprintf('%s/%s', refPath, [srcName '_DS.tif']);
    if exist(srcDsFileName, 'file')
        delete(srcDsFileName);
    end
    fprintf('Writing %s\n', srcDsFileName);
    dos(sprintf('gdalwarp -srcnodata "0" -dstnodata "0" -r cubicspline -tap -tr %d %d "%s" "%s"', ...        
        abs(refR.DeltaX), abs(refR.DeltaY), srcFileName, srcDsFileName));

    [srcDsIm srcDsR] = geotiffread(srcDsFileName);
%     srcFileName = srcDsFileName;
else
    exit;
end

if (~all(size(refIm) == size(srcDsIm)))
    %crop out relevant section of ref image
    [x, y] = srcDsR.intrinsicToWorld(1,1);
    [ulJ, ulI] = refR.worldToIntrinsic(x, y);

    [x, y] = srcDsR.intrinsicToWorld(size(srcDsIm, 2), size(srcDsIm, 1));
    [brJ, brI] = refR.worldToIntrinsic(x, y);

    refSubIm = refIm(ulI:brI, ulJ:brJ, :);
%     refR = 
else
    refSubIm = refIm;
end

figure;
subplot(1,3,1)
imagesc(srcDsIm(:,:,[1:3])*16)
subplot(1,3,2)
imagesc(uint16(refSubIm(:,:,[1:3])*16))
subplot(1,3,3)
imshow(uint16(refIm(:,:,[1:3])*16))
hold on;
rectangle('Position', [ulJ ulI  size(srcDsIm, 2) size(srcDsIm, 1)], 'EdgeColor', 'r')

% refIm = modisSubIm;
% srcDsIm = ngiDsIm;
if true
    numBands = size(srcDsIm, 3);
    for i = 1:size(srcDsIm, 1)
        rowIdx = max(1,i-(winSize-1)/2):min(size(refSubIm, 1), i+(winSize-1)/2);
        for j = 1:size(srcDsIm, 2)
            colIdx = max(1,j-(winSize-1)/2):min(size(refSubIm, 2), j+(winSize-1)/2);
            refWin = reshape(refSubIm(rowIdx, colIdx, :), [], numBands);
            srcWin = reshape(srcDsIm(rowIdx, colIdx, :), [], numBands);
            tmp = double(refWin)./double(srcWin);
            tmp(tmp==inf)=nan;
            gainIm(i, j, :) = nanmean(tmp, 1);
        end
    end

%average gain over bands
% gainImMean = nanmean(gainIm, 3);
% gainImMean = gainImMean./mean(gainImMean(:));
% for b = 1:numBands
%     band = gainIm(:, :, b);
%     bandMean(b) = nanmean(band(:));
%     gainIm(:, :, b) = bandMean(b) * gainImMean;
% end
else  %find spatially varying gain from intensity ratio
    numBands = size(srcDsIm, 3);
    for i = 1:size(srcDsIm, 1)
        rowIdx = max(1,i-(winSize-1)/2):min(size(refSubIm, 1), i+(winSize-1)/2);
        for j = 1:size(srcDsIm, 2)
            colIdx = max(1,j-(winSize-1)/2):min(size(refSubIm, 2), j+(winSize-1)/2);
            refWin = reshape(refSubIm(rowIdx, colIdx, :), [], numBands);
            srcWin = reshape(srcDsIm(rowIdx, colIdx, :), [], numBands);
            %prevent nodata contributing to mean vals
            refWin(srcWin == 0) = nan;
            srcWin(srcWin == 0) = nan;
            tmp = nanmean(double(refWin), 2)./nanmean(double(srcWin), 2);
            tmp(tmp == inf) = nan;
%             tmp(tmp == 0) = nan;
            gainIm(i, j, :) = nanmean(tmp, 1);
        end
    end

     gainImMean = gainIm;
    gainImMean = gainImMean./nanmean(gainImMean(:));
    for b = 1:numBands
        srcBand = srcDsIm(:, :, b);
        refBand = refSubIm(:, :, b);
        %prevent nodata contributing to mean vals        
        refBand(srcBand == 0) = nan;
        srcBand(srcBand == 0) = nan;
        srcBandMean(b) = nanmean(srcBand(:));
        refBandMean(b) = nanmean(refBand(:));
        gainIm(:, :, b) = gainImMean * refBandMean(b) / srcBandMean(b);
    end
end

gainDsFileName = sprintf('%s/%s', refPath, [srcName '_DS_GAIN.tif']);
if exist(gainDsFileName, 'file')
     delete(gainDsFileName);
end
%TO DO: no data?
uint16ScaleFactor = (2^12-1);
srcDsInfo = geotiffinfo(srcDsFileName);
tmp = srcDsInfo.GeoTIFFTags.GeoKeyDirectoryTag;
tmp.ProjFalseEastingGeoKey = 21;
fprintf('Writing %s\n', gainDsFileName);
geotiffwrite(gainDsFileName, uint16(gainIm*uint16ScaleFactor), srcDsInfo.RefMatrix, 'GeoKeyDirectoryTag', ...
    tmp);
% i = geotiffinfo(gainDsFileName);i.ProjParm

%fix the mapping toolbox's bug with long origin
%the PROJCS... string below was obtained with 'gdalsrsinfo -p -o wkt_simple ...
fprintf('Updating %s\n', gainDsFileName);
dos(sprintf('gdal_edit.py -a_nodata 0 -a_srs "%s" "%s"', ...
    'PROJCS["unnamed", GEOGCS["unknown", DATUM["unknown", SPHEROID["WGS84",6378137,298.257223563]], PRIMEM["Greenwich",0], UNIT["degree",0.0174532925199433]], PROJECTION["Transverse_Mercator"], PARAMETER["latitude_of_origin",0], PARAMETER["central_meridian",21], PARAMETER["scale_factor",1], PARAMETER["false_easting",0], PARAMETER["false_northing",0], UNIT["metre",1]]', ...
    gainDsFileName));

% gainIm_ = gainIm;
% gainIm_(gainIm_ == inf) = 0;
% 
% gainIm_ = gainIm_/max(gainIm_(:));

[srcIm srcR] = geotiffread(srcFileName);

%upsample gains and apply to orig image
if (true)
    %use GDAL
%     refPath = fileparts(refFileName);
%     [srcDsPath srcDsName] = fileparts(srcDsFileName);
%     gainDsFileName = sprintf('%s/%s', refPath, [srcDsName '_GAIN.tif']);
%     if exist(gainDsFileName, 'file')
%          delete(gainDsFileName);
%     end
%     %TO DO: no data?
%     uint16ScaleFactor = (2^12-1);
%     srcDsInfo = geotiffinfo(srcDsFileName);
%     geotiffwrite(gainDsFileName, uint16(gainIm*uint16ScaleFactor), srcDsInfo.SpatialRef, 'GeoKeyDirectoryTag', ...
%         srcInfo.GeoTIFFTags.GeoKeyDirectoryTag);

    [srcPath srcName] = fileparts(srcFileName);
    gainUsFileName = sprintf('%s/%s', refPath, [srcName '_US_GAIN.tif']);
    if exist(gainUsFileName, 'file')
         delete(gainUsFileName);
    end
%     dos(sprintf('gdalwarp -multi -r lanczos -tap -tr %d %d "%s" "%s"', ...        
%         abs(srcInfo.SpatialRef.DeltaX), abs(srcInfo.SpatialRef.DeltaY), gainDsFileName, gainUsFileName));
    dos(sprintf('gdalwarp -srcnodata "0" -dstnodata "0" -wm 1024 -multi -r cubicspline -ts %d %d "%s" "%s"', ...        
        abs(srcInfo.SpatialRef.RasterSize(2)), abs(srcInfo.SpatialRef.RasterSize(1)), gainDsFileName, gainUsFileName));
    [gainImUs gainR] = geotiffread(gainUsFileName);
    
    gainImUs = (single(gainImUs)/single(uint16ScaleFactor));
    calibIm = uint16(gainImUs.*single(srcIm));
    gainIm = gainImUs;
    
else
    %use Matlab
%     [ngiIm ngiR] = geotiffread(ngiDsFileName2);
%     srcNgiIm = ngiIm;
    gainImUs = zeros(size(srcIm));
    calibIm = srcIm;
    for i = 1:4
        fprintf('.');
        band = single(srcIm(:,:,i));
        gainImUs(:,:,i) = imresize(gainIm(:,:,i), size(band), 'lanczos2'); %'lanczos2'
        calibBand = uint16(single(band).*gainImUs(:,:,i));
        calibIm(:,:,i) = calibBand;
    end
    fprintf('\n');
    gainIm = gainImUs;
    
end
if (true)
    refPath = fileparts(refFileName);
    [srcPath srcName] = fileparts(srcFileName);
    calibFileName = sprintf('%s/%s', refPath, [srcName '_XCALIB.tif']);
    if exist(calibFileName, 'file')
         delete(calibFileName);
    end
    fprintf('Writing %s\n', calibFileName);        
    geotiffwrite(calibFileName, calibIm, srcR, 'GeoKeyDirectoryTag', ...
        srcInfo.GeoTIFFTags.GeoKeyDirectoryTag);
    dos(sprintf('gdal_edit.py -a_nodata 0 -a_srs "%s" "%s"', ...
        'PROJCS["unnamed", GEOGCS["unknown", DATUM["unknown", SPHEROID["WGS84",6378137,298.257223563]], PRIMEM["Greenwich",0], UNIT["degree",0.0174532925199433]], PROJECTION["Transverse_Mercator"], PARAMETER["latitude_of_origin",0], PARAMETER["central_meridian",21], PARAMETER["scale_factor",1], PARAMETER["false_easting",0], PARAMETER["false_northing",0], UNIT["metre",1]]', ...
        calibFileName));

end

[refStretchIm_ refStretchLim] = StretchImage(refIm);
refStretchIm = StretchImage(refSubIm, 'stretchLim', refStretchLim);
srcStretchIm = StretchImage(srcIm, 'stretchLim', refStretchLim);
gainIm_ = gainIm./max(gainIm(:));
calibStretchIm = StretchImage(calibIm, 'stretchLim', refStretchLim);

Plot(refStretchIm, gainIm_, srcStretchIm, calibStretchIm, 1:3);
Plot(refStretchIm, gainIm_, srcStretchIm, calibStretchIm, [4 1 2]);

res.refIm  = refIm;
res.refR  = refR;
res.refSubIm  = refSubIm;
res.refSubR = srcDsR;
res.srcIm = srcIm;
res.calibIm = calibIm;
res.calibR = srcR;
res.gainIm = gainIm;
res.refStretchIm = refStretchIm;
res.srcStretchIm = srcStretchIm;
res.calibStretchIm = calibStretchIm;
end

function Plot(refStretchIm, gainIm, srcStretchIm, calibStretchIm, bands)
figure;
subplot(2,2,1);
imagesc((refStretchIm(:,:,bands)));
title('Ref')
axis image off
subplot(2,2,2);
imagesc(gainIm(:, :, bands));
title('Gain')
axis image off
h1 = subplot(2,2,3);
imagesc((srcStretchIm(:,:,bands)));
title('Source')
axis image off
h2 = subplot(2,2,4);
imagesc(calibStretchIm(:,:,bands));
title('Calib')
axis image off

linkaxes([h1 h2], 'xy')
end

