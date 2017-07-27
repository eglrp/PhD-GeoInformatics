%% ------------------------------------------------------------------------------------ 
% Make MODIS tiff from hdf in NGI projection

%Naming convention
%-----------------
% MOD09A1.A2006001.h08v05.005.2006012234657.hdf indicates:
% 
%     MOD09A1 - Product Short Name
%     .A2006001 - Julian Date of Acquisition (A-YYYYDDD)
%     .h08v05 - Tile Identifier (horizontalXXverticalYY)
%     .005 - Collection Version
%     .2006012234567 - Julian Date of Production (YYYYDDDHHMMSS)
%     .hdf - Data Format (HDF-EOS)
%
% BANDS 
%-----------------
% 1 	620–670 	  	Absolute Land Cover Transformation, Vegetation Chlorophyll
% 2 	841–876 	  	Cloud Amount, Vegetation Land Cover Transformation
% 3 	459–479 	  	Soil/Vegetation Differences
% 4 	545–565 	  	Green Vegetation
% 5 	1230–1250 	  	Leaf/Canopy Differences
% 6 	1628–1652 	  	Snow/Cloud Differences
% 7 	2105–2155 	  	Cloud Properties, Land Properties
%
% Make RGBN image from bands 1,4,3,2
%DMC bands
% Blue: 400-580 nm n Green: 500-650 nm n Red: 590-675 nm n Near infrared: 675-850 nm n Near infrared alternate: 740-850 nm
%
% MISR 
% The center wavelength of each of these bands is 446, 558, 672, and 867 nanometers respectively.

%MODIS bands
%1 	620 - 670
%4 	545 - 565
%3 	459 - 479
%2 is 841 - 876

%SPOT bands
% 	Panchromatic, 2.5/5m, 0.48 - 0.71 µm
% 	B1 : green, 10m, 0.50 - 0.59 µm
% 	B2 : red, 10m, 0.61 - 0.68 µm
% 	B3 : near infrared, 10m, 0.78 - 0.89 µm
% 	B4 : mid infrared (MIR), 20m, 1.58 - 1.75 µm

close all;clear all;
modisFilename = 'F:\MSc GeoInformatics\Data\MODIS\MCD43A4.A2010025.h19v12.005.2010043064233.hdf';
% modisFilename = 'F:\MSc GeoInformatics\Data\MODIS\MOD09A1.A2010033.h19v12.005.2010043033652.hdf';
% ngiFileName = 'G:\MSc GeoInformatics\Data\NGI\Rectified\RGB\3320AA\3320AA_01_2010_312_RGB_RECT.tif';
pciFileName = 'F:\MSc GeoInformatics\Data\NGI\Cross Calibration\SmallTest\3321B_3172_01_0018_RGBN_DS_2_5.tif';

% dos(sprintf('gdalinfo -proj4 "%s"', ngiFileName)) %'+proj=tmerc +lat_0=0 +lon_0=21 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs ' 
dos(sprintf('gdalinfo -proj4 "%s"', modisFilename)) %+proj=tmerc +lat_0=0 +lon_0=21 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs '
dos(sprintf('gdalinfo -proj4 "%s"', pciFileName)) %+proj=tmerc +lat_0=0 +lon_0=21 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs '

%Warp modis bands to tiff files in pci/ngi projection
bands = [1 4 3 2];
[d n] = fileparts(modisFilename);
[p d] = fileparts(pciFileName);
for i = 1:length(bands)
    outFileName{i} = sprintf('%s/%s.Lo21.Band%d.tif', p, n, bands(i));
    if exist(outFileName{i}, 'file')
        delete(outFileName{i});
    end
    
    dos(sprintf('gdalwarp -r lanczos -tr 500 500 -tap -t_srs "%s" HDF4_EOS:EOS_GRID:"%s":MOD_Grid_BRDF:Nadir_Reflectance_Band%d "%s"', ...        
        '+proj=tmerc +lat_0=0 +lon_0=21 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs', modisFilename, ...
        bands(i), outFileName{i}));
%     dos(sprintf('gdalwarp -r lanczos -tr 500 500 -tap -t_srs "%s" HDF4_EOS:EOS_GRID:"%s":MOD_Grid_500m_Surface_Reflectance:sur_refl_b0%d "%s"', ...        
%         '+proj=tmerc +lat_0=0 +lon_0=23 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs', modisFilename, ...
%         bands(i), outFileName{i}));
end

%Read tiff modis band images into one image
rgbnOutFileName = sprintf('%s/%s.Lo21.RGBN.tif', p, n);

if exist(rgbnOutFileName, 'file')
    delete(rgbnOutFileName);
end

if false
    %NB combining bands with mapping toolbox messes up projection info
    %somehow
    for i = 1:length(bands)
        [modisIm(:,:,(i)) r{i}] = geotiffread(outFileName{i});
    end

    pciInfo = geotiffinfo(pciFileName)
    modisInfo = geotiffinfo(outFileName{1})


    geotiffwrite(rgbnOutFileName, modisIm, r{end}, 'GeoKeyDirectoryTag', ...
        pciInfo.GeoTIFFTags.GeoKeyDirectoryTag);
else
% gdal_merge.py -seperate red.tif nir.tif mir.tif -o stack.tif
    bandFileStr = '';
    for i = 1:length(bands)
        bandFileStr = sprintf('%s "%s"', bandFileStr, outFileName{i});
    end
    dos(sprintf('gdal_merge.py -seperate %s -o "%s"', bandFileStr, rgbnOutFileName));

end

modisRgbnInfo = geotiffinfo(rgbnOutFileName)

[testIm testR] = geotiffread(rgbnOutFileName);

clear im;
for i = 1:length(bands)
    data = testIm(:,:,i);
    nonodata = data(data~=(2^15-1));
    nonodata = double(nonodata)./(2^15-1);
    data(data==(2^15-1))=0;
    data = double(data)./(2^15-1);
    
    im(:,:,i) = imadjust(data, stretchlim(nonodata,[0 .99]), [0 1], 1);    
end

rgbIm = im(:,:,1:3);
cirIm = im(:,:,[4 1 2]);

figure;
subplot(1,2,1)
imagesc(rgbIm)
subplot(1,2,2)
imagesc(cirIm)

%--------------------------------------------------------------------------
%%
clear all close all

modisFilename = 'G:\MSc GeoInformatics\Data\MODIS\MCD43A4.A2010025.h19v12.005.201004306.RGBN.tif';
ngiFileName = 'H:\MSc GeoInformatics\Data\NGI\My Rectified\3321D_2010_319\o3321D_2010_319_01_0004_RGBIR.tif';

%downsample the NGI image with GDAL to MODIS res
[p n e] = fileparts(modisFilename);
ngiDsFileName = sprintf('%s/%s', p, 'ngiDs250.tif');
delete(ngiDsFileName);
dos(sprintf('gdalwarp -r bilinear -tap -tr 250 250 "%s" "%s"', ...        
    ngiFileName, ngiDsFileName));

%downsample the NGI image with GDAL to workable high res
ngiDsFileName2 = sprintf('%s/%s', p, 'ngiDs5.tif');
delete(ngiDsFileName2);
dos(sprintf('gdalwarp -r bilinear -tap -tr 2.5 2.5 "%s" "%s"', ...        
    ngiFileName, ngiDsFileName2));

%read images
[modisIm modisR] = geotiffread(modisFilename);
% [ngiIm ngiR] = geotiffread(ngiFileName);

[ngiDsIm ngiDsR] = geotiffread(ngiDsFileName);


%crop out relevant section of MODIS image
[x, y] = ngiDsR.intrinsicToWorld(1,1);
[ulJ, ulI] = modisR.worldToIntrinsic(x, y);

[x, y] = ngiDsR.intrinsicToWorld(size(ngiDsIm, 2), size(ngiDsIm, 1));
[brJ, brI] = modisR.worldToIntrinsic(x, y);

modisSubIm = modisIm(ulI:brI, ulJ:brJ, :);

figure;
subplot(1,3,1)
imagesc(ngiDsIm(:,:,[1:3])*16)
subplot(1,3,2)
imagesc(uint16(modisSubIm(:,:,[1:3])*16))

figure
subplot(1,2,1)
% imshow(uint16(ngiIm(:,:,[1:3])*16))
subplot(1,2,2)
imshow(uint16(modisIm(:,:,[1:3])*16))
hold on;
rectangle('Position', [ulJ ulI  size(ngiDsIm, 2) size(ngiDsIm, 1)], 'EdgeColor', 'r')

%%
% fun = @(x) median(x(:));
% B = nlfilter(modisIm,[3 3 4],fun);

refIm = modisSubIm;
calibIm = ngiDsIm;
winSize = 3;
numBands = 4;
for i = 1:size(ngiDsIm, 1)
    rowIdx = max(1,i-(winSize-1)/2):min(size(refIm, 1), i+(winSize-1)/2);
    for j = 1:size(ngiDsIm, 2)
        colIdx = max(1,j-(winSize-1)/2):min(size(refIm, 2), j+(winSize-1)/2);
        refWin = reshape(refIm(rowIdx, colIdx, :), [], numBands);
        calibWin = reshape(calibIm(rowIdx, colIdx, :), [], numBands);
        tmp = double(refWin)./double(calibWin);
        tmp(tmp==inf)=nan;
%         tmp = tmp(~isnan(tmp));
        gainIm(i, j, :) = nanmean(tmp, 1);
    end
end

gainIm_ = gainIm;
gainIm_(gainIm_ == inf) = 0;

gainIm_ = gainIm_/max(gainIm_(:));

figure;
subplot(1,3,1)
imagesc(uint16(16*refIm(:,:,1:3)))
subplot(1,3,2)
imagesc(uint16(calibIm(:,:,1:3)*16))
subplot(1,3,3)
imagesc(gainIm_(:,:,1:3))


%%
%upsample gaimIm and apply

% gainIm = imresize(gainIm, [size(ngiIm, 1) size(ngiIm, 2)], 'bicubic');
% imagesc(gainIm(1:10:end,1:10:end,1))
% calibNgiIm = ngiIm;
[ngiIm ngiR] = geotiffread(ngiDsFileName2);
calibNgiIm = ngiIm;
gainImUs = zeros(size(ngiIm));
for i = 1:4
    fprintf('.');
    band = single(ngiIm(:,:,i));
    gainImUs(:,:,i) = imresize(gainIm(:,:,i), size(band), 'lanczos2'); %'lanczos2'
    calibBand = uint16(single(band).*gainImUs(:,:,i));
    calibNgiIm(:,:,i) = calibBand;
end

[p n e] = fileparts(modisFilename);
ngiCalibFileName = sprintf('%s/%s', p, 'ngiCalib.tif');

ngiInfo = geotiffinfo(ngiFileName);
geotiffwrite(ngiCalibFileName, calibNgiIm, ngiR, 'GeoKeyDirectoryTag', ...
    ngiInfo.GeoTIFFTags.GeoKeyDirectoryTag);

figure;
subplot(1,3,1)
imshow(uint16(modisSubIm(:,:,1:3)*16))
h1=subplot(1,3,2);
imshow(ngiIm(1:10:end,1:10:end,1:3)*16)
h2=subplot(1,3,3);
imshow(calibNgiIm(1:10:end,1:10:end,1:3)*16)
linkaxes([h1 h2], 'xy');

figure;
imshow(gainImUs(:,:,1:3))


%%
%DO the same as the above but using backend function
clear all; 
close all

% modisFilename = G:\MSc GeoInformatics\Data\NGI\Cross Calibration\SmallTest\MCD43A4.A2010025.h19v12.005.201004306.RGBN.tif';
modisFilename = 'G:\MSc GeoInformatics\Data\NGI\Cross Calibration\SmallTest/MOD09A1.A2010033.h19v12.005.2010043033652.RGBN.tif';
% modisFilename = 'G:\MSc GeoInformatics\Data\NGI\Cross Calibration\SmallTest\MOD09A1.A2010033.h19v12.005.201004303.RGBN.tif';
% ngiFileName = 'H:\MSc GeoInformatics\Data\NGI\My Rectified\3321D_2010_319\o3321D_2010_319_01_0004_RGBIR.tif';
ngiFileName = 'H:\MSc GeoInformatics\Data\NGI\My Rectified\3321D_2010_319\o3321D_2010_319_01_0022_RGBIR.tif';
% ngiFileName = 'H:\MSc GeoInformatics\Data\NGI\My Rectified\3321D_2010_319\o3321D_2010_319_02_0051_RGBIR.tif';

[pModis] = fileparts(modisFilename);
%downsample the NGI image with GDAL to more workable high res
[pNgi nNgi] = fileparts(ngiFileName);
ngiDsFileName = sprintf('%s/%s', pModis, [nNgi '_Ds_2pt5.tif']);
delete(ngiDsFileName);
dos(sprintf('gdalwarp -r bilinear -tap -tr 2.5 2.5 "%s" "%s"', ...        
    ngiFileName, ngiDsFileName));

gainIm2 = CrossCalib(modisFilename, ngiDsFileName);

%% As Above but with already DS files from SunGis08
% clear all; 
close all;

modisFilename = 'F:\MSc GeoInformatics\Data\NGI\Cross Calibration\SmallTest\MCD43A4.A2010025.h19v12.005.2010043064233.RGBN.tif';
% modisFilename = 'F:\MSc GeoInformatics\Data\NGI\Cross Calibration\SmallTest\MOD09A1.A2010033.h19v12.005.2010043033652.RGBN.tif';

% ngiDsFileName = {'F:\MSc GeoInformatics\Data\NGI\Cross Calibration\SmallTest\3321B_3172_01_0018_RGBN_DS_2_5.tif';...
%     'F:\MSc GeoInformatics\Data\NGI\Cross Calibration\SmallTest\3321B_3172_02_0053_RGBN_DS_2_5.tif';};
% ngiDsFileName = {'F:\MSc GeoInformatics\Data\NGI\Cross Calibration\SmallTest\3321D_319_01_0022_RGBN_DS_2_5.tif';...
%     'F:\MSc GeoInformatics\Data\NGI\Cross Calibration\SmallTest\3321D_319_02_0051_RGBN_DS_2_5.tif';...
%     'F:\MSc GeoInformatics\Data\NGI\Cross
%     Calibration\SmallTest\3321B_3172_12_0416_RGBN_DS_2_5.tif'...
% };
ngiDsFileName = {'F:\MSc GeoInformatics\Data\NGI\Cross Calibration\SmallTest\3321b_3172_01_0025_rgbn_INIT.tif';...
    'F:\MSc GeoInformatics\Data\NGI\Cross Calibration\SmallTest\3321b_3172_01_0026_rgbn_INIT.tif'...
    };
ngiDsFileName = {'F:\MSc GeoInformatics\Data\NGI\Cross Calibration\SmallTest\3321b_3172_05_0159_rgbn_INIT.tif';...
    'F:\MSc GeoInformatics\Data\NGI\Cross Calibration\SmallTest\3321b_3172_05_0160_rgbn_INIT.tif'...
    };

% figure;
for i = 1:length(ngiDsFileName)
    %[refIm{i} srcIm{i} calibIm{i} refSpatialRef{i} calibSpatialRef{i} gainIm{i}] = CrossCalib(modisFilename, ngiDsFileName{i});
    res(i) = CrossCalib(modisFilename, ngiDsFileName{i});
end

%%
f1 = figure;
f2 = figure;
calibStretchLim = [];
srcStretchLim = [];
for i = 1:length(ngiDsFileName)
%     if (isempty(calibStretchLim))
%         [calibStretchIm{i} calibStretchLim] = StretchImage(res(i).calibIm);
%         [refStretchIm{i}] = StretchImage(uint16(res(i).refIm), 'stretchLim', calibStretchLim);
%     else
%         calibStretchIm{i} = StretchImage(res(i).calibIm, 'stretchLim', calibStretchLim);
%         [refStretchIm{i}] = StretchImage(uint16(res(i).refIm), 'stretchLim', calibStretchLim);
%     end
%     if (isempty(srcStretchLim))
%         [srcStretchIm{i} srcStretchLim] = StretchImage(res(i).srcIm);
%     else
%         srcStretchIm{i} = StretchImage(res(i).srcIm, 'stretchLim', srcStretchLim);
%     end

    figure(f1)
    subplot(1,3,1)
%     si = geotiffinfo(modisFilename);
%     modisSpatialRef = spatialRef{i};
    mapshow(res(i).refStretchIm(:,:,1:3), res(i).refSubR);
    axis tight 
    axis off
    hold on;
    subplot(1,3,2)
    mapshow(res(i).srcStretchIm(:,:,1:3), res(i).calibR);
    axis tight 
    axis off
    hold on;
    subplot(1,3,3)
    mapshow(res(i).calibStretchIm(:,:,1:3), res(i).calibR);
    axis tight 
    axis off
    hold on;

    figure(f2)
    subplot(1,3,1)
    mapshow(res(i).refStretchIm(:,:,[4 1 2]), res(i).refSubR);
    axis tight 
    axis off
    hold on;
    subplot(1,3,2)
    mapshow(res(i).srcStretchIm(:,:,[4 1 2]), res(i).calibR);
    axis tight 
    axis off
    hold on;
    subplot(1,3,3)
    mapshow(res(i).calibStretchIm(:,:,[4 1 2]), res(i).calibR);
    axis tight 
    axis off
    hold on;
end

%%
%Colour balance source images to match calib (ref) images

[colourBalImage whiteBal] = ColourBalImage(refIm, srcIm, varargin)

%%
% Visualise gain images
clear all; 
close all;

% modisFilename = 'F:\MSc GeoInformatics\Data\NGI\Cross Calibration\SmallTest\MCD43A4.A2010025.h19v12.005.201004306.RGBN.tif';
modisFilename = 'F:\MSc GeoInformatics\Data\NGI\Cross Calibration\SmallTest\MOD09A1.A2010033.h19v12.005.2010043033652.RGBN.tif';

% gainDsFileName = {'F:\MSc GeoInformatics\Data\NGI\Cross Calibration\SmallTest\3321B_3172_01_0018_RGBN_DS_2_5_DS_GAIN.tif';...
%     'F:\MSc GeoInformatics\Data\NGI\Cross Calibration\SmallTest\3321B_3172_02_0053_RGBN_DS_2_5_DS_GAIN.tif';};
gainDsFileName = {'F:\MSc GeoInformatics\Data\NGI\Cross Calibration\SmallTest\3321D_319_01_0022_RGBN_DS_2_5_DS_GAIN.tif';...
    'F:\MSc GeoInformatics\Data\NGI\Cross Calibration\SmallTest\3321D_319_02_0051_RGBN_DS_2_5_DS_GAIN.tif';...
    'F:\MSc GeoInformatics\Data\NGI\Cross Calibration\SmallTest\3321B_3172_12_0416_RGBN_DS_2_5_DS_GAIN.tif'};

figure;
for i = 1:length(gainDsFileName)
    [gainIm{i} r] = geotiffread(gainDsFileName{i});
    subplot(1, length(gainDsFileName), i)
    imagesc(5*gainIm{i}(:,:,1:3));
end

for i = 1:length(gainDsFileName)
    figure;
    for j = 1:4
        band = double(gainIm{i}(:,:,j));
        band = band./max(band(:));
        subplot(2, 2, j)
        imagesc(band);
    end
end

PlotMultibandImages(gainIm)
PlotMultibandImage(gainIm{1})
%%
% Visualise gain images (compare c++ & matlab)
clear all; 
close all;

% modisFilename = 'F:\MSc GeoInformatics\Data\NGI\Cross Calibration\SmallTest\MCD43A4.A2010025.h19v12.005.201004306.RGBN.tif';
modisFilename = 'F:\MSc GeoInformatics\Data\NGI\Cross Calibration\SmallTest\MOD09A1.A2010033.h19v12.005.2010043033652.RGBN.tif';

% gainDsFileName = {'F:\MSc GeoInformatics\Data\NGI\Cross Calibration\SmallTest\3321B_3172_01_0018_RGBN_DS_2_5_DS_GAIN.tif';...
%     'F:\MSc GeoInformatics\Data\NGI\Cross Calibration\SmallTest\3321B_3172_02_0053_RGBN_DS_2_5_DS_GAIN.tif';};
gainDsFileName = {'F:\MSc GeoInformatics\Data\NGI\Cross Calibration\SmallTest\3321b_3172_05_0160_rgbn_INIT_DS_GAIN.tif';...
    'F:\MSc GeoInformatics\Data\NGI\Cross Calibration\SmallTest\3321b_3172_05_0160_rgbn_INIT_INIT_DS_GAIN.tif';...
    'F:\MSc GeoInformatics\Data\NGI\Cross Calibration\SmallTest\3321b_3172_05_0159_rgbn_INIT_DS_GAIN.tif';...
    'F:\MSc GeoInformatics\Data\NGI\Cross Calibration\SmallTest\3321b_3172_05_0159_rgbn_INIT_INIT_DS_GAIN.tif'...
    };

figure;
for i = 1:length(gainDsFileName)
    [gainIm{i} r] = geotiffread(gainDsFileName{i});
    subplot(2, 2, i)
    if (class(gainIm{i}) == 'uint16')
        gainIm{i} = double(gainIm{i})/(2^12 - 1);
    end
    imagesc(gainIm{i}(:,:,1:3)/5);
end

%NOTES
%---------------------------------------------------------------------
%- Looking at the mosaic of a job on SunGis08, there are numerous
%discontinuities between images.
%- Closer inspection shows that it is on the edges of the overlap region
%that this occurs, deeper into the overlap region, the join is smoother.
%- What is the gdal_warp behaviour when downsampling and there are
%partially covered pixels on the edges?  Are these set to a value or to
%nodata.  If they were set to nodata, we would expect to lose data when
%upsampling again as the partial pixels would have been dropped.  They are
%set to a value seemingly if they are greater than half the area of a full
%pixel otherwise set to nodata. We need to exclude them somehow.
%- It shouldn't be necessary to use the same upsampling as downsampling
%technique.  But we should make sure the downsampling is averaging the
%smaller pixels to make the bigger ones, not just taking the center calue
%or something.  So perhaps use 'average'

%%
%Compare matlab and c++ gain ims
close all;
gainDsFileName = {...
    'F:\MSc GeoInformatics\Data\NGI\Cross Calibration\SmallTest\3321D_319_01_0022_RGBN_DS_2_5_DS_GAIN.tif';...
    'F:\MSc GeoInformatics\Data\NGI\Cross Calibration\SmallTest\3321D_319_01_0022_RGBN_DS_2_5_DS_GAIN2.tif';...
    };

uint16ScaleFactor = (2^12-1);
figure;

[gainIm r] = geotiffread(gainDsFileName{1});
gainIm1 = double(gainIm(:, :, 1:3))/(5*uint16ScaleFactor);
subplot(1,3,1)
imagesc(gainIm1);

[gainIm r] = geotiffread(gainDsFileName{2});
gainIm2 = double(gainIm(:, :, 1:3))/(5);
subplot(1,3,2)
imagesc(gainIm2);

subplot(1,3,3)
imagesc(abs(gainIm2 - gainIm1)*5);

% gainIm(gainIm==0) = nan;
% PlotMultibandImage(gainIm);

%%
fn = 'F:\MSc GeoInformatics\Data\NGI\Cross Calibration\SmallTest\3321D_319_02_0051_RGBN_DS_2_5_INIT_US_GAIN.tif';
[im r] = geotiffread(fn);

figure
imshow(im(:,:,1:3)/5)

fn = 'F:\MSc GeoInformatics\Data\NGI\Cross Calibration\SmallTest\3321D_319_02_0051_RGBN_DS_2_5_INIT_XCALIB.tif';
[im r] = geotiffread(fn);

figure
imshow(im(:,:,1:3)*16)

%%
%Notes (on gain ims)
%-------------------
% - The gain between bands is fairly similar and may point to the
% possibility of averaging over bands.
% - The gain between images (in same job) is vaguely similar but shows
% correlation with patterns in the images themselves which is problematic.
% - There is a case for smoothing gains and or averaging over bands and or
% more than one gain image
% - If averaging over more than one gain image, it would be a lot easier to
% do before rectification.  Not sure if this makes sense though - part of
% the effect is due to sun/camera pos and part is due to land cover.  Land
% cover varies over images so probably this doesn't make sense.  If we are
% only correcting from NGI im to matching area in MODIS, it should be done
% after rectification as this is the only way we can match the areas.
% - A test of averaging gains over bands looks ok.  Bit of a difference vis
% between 2 images but looks pretty good.
% - A test of gains from intensity looks ok.  Worse differences between
% adjacent images than above but not bad.  Note (sum(refBands)/sum(srcBands))
% does not = sum(refBands/srcBands)


%%
% TO DO
%-------
% - Make small before and after calib mosaic - hopefully across some
% boundary that shows it is helping (would be nice to show MODIS mosaic
% that shows recognisable feature)
%- Download the above ims from SunGis08
%- See what the effect of the different calib params are: 
%  - interp method,
%  - modis res (pan sharpen to 250m?), 
%  - window size
%  - MODIS im (8 day/ 16 day /pan sharpened)
%  - gain + offset
%  - intensity gain only
%- Compare surfaces of calib gains - should show some kind of pattern? eg
%BRDF
%- How much work would it be to write a GDAL exe to find gain, upsample and
%apply?  Then we could run that on the server.
%- Do we have the correct MODIS image? Date, bands
%- Can using another MODIS image eg 8 day help?
%- Can MODIS be pan sharpened to 250m using bands 1&2
%- Try matching only 1 band betw NGI and MODIS i.e. the one that is
%spectrally the most similar and then use this as a general intensity
%adjustment.  It it happens to be MODIS band 1 or 2 then use one of the
%250m options.  The closest looks like the red.
%- Smooth the high res (have some smoothness constraint) rather than
%interpolate
%- Error measures like downsampling calibrated im back to MODIS res then
%subtracting from ref MODIS im


%New TO DO
%----------
%- Look at gain images for some kind of pattern (within 1 job and close 
%together in time).  Bear in mind that exposure differences will disturb
%the pattern.
%- Should radiometric calib be done before or after geometric calib?
%- Maybe gain images should averaged over close times.  Like another moving
%window approach but in time (& space).  Would need to deal with different
%imsizes due to georectification.
%- Compare gains in different bands - do they match?  Can we use intensity
%only to find the gains?

