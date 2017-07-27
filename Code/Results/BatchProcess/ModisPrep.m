%NB Note: WAMIS images are not calibrated to be used for visuals only - I
%should use MOSIS NBAR 8 day composite or MOD09, available from NASA
%website

%make 4 band MODIS image with bands roughly corresp to NGI
%DMC bands
% Blue: 400-580 nm n Green: 500-650 nm n Red: 590-675 nm n Near infrared: 675-850 nm n Near infrared alternate: 740-850 nm
%
%MODIS bands
%1 	620 - 670
%4 	545 - 565
%3 	459 - 479
%2 is 841 - 876

%Note that the MODIS bands match the DMC ones pretty poorly which makes
%this calibration suspect purely from this spectral perspective

fileNames = {'F:\MSc GeoInformatics\Data\MODIS\t1.2010040.0755.utm35j-2.SouthAfrica.143.231m.paa.tif', ...
    'F:\MSc GeoInformatics\Data\MODIS\t1.2010040.0755.utm35j-2.SouthAfrica.721.231m.paa.tif'};

for i = 1:length(fileNames)
    t{i} = Tiff(fileNames{i});
    im{i} = t{i}.read();
%     t{i}.close();
end

tout = Tiff('F:\MSc GeoInformatics\Data\MODIS\t1.2010040.0755.utm35j-2.SouthAfrica.1432.231m.paa.tif','w');

tagId = Tiff.TagID;
fn = fieldnames(tagId);
tagStruct = tagId;
for i = 1:length(fn)
    try
        val = t{1}.getTag(fn{i});
        tout.setTag(fn{i}, val);
        fprintf('Setting tag "%s" to %d\n', fn{i}, val);        
    catch
    end
end

imout = im{1};
imout = cat(3, imout, im{2}(:,:,2));
tout.setTag('Photometric', 1);
% tout.setTag('ExtraSamples', 0);
tout.setTag('BitsPerSample', 8);
tout.setTag('SamplesPerPixel', 4);
tout.setTag('MaxSampleValue', 255);
% tout.setTag('MaxSampleValue', 4095);
% tout.setTag('NoDataValue', 4095);
tout.setTag('Compression', Tiff.Compression.LZW); %DEFLATE IS THE PNG ALGORITHM
% tout.rewriteDirectory();
tout.write(imout);
% tout.setTag('BitsPerSample', 12);
% tout.rewriteDirectory();
tout.close();
t{1}.close();
t{2}.close();

%%
% Look for other files that could be breaking the mosaic in RGB

baseDir = 'F:\MSc GeoInformatics\Data\NGI\Rectified\RGB\';
dirList = dir([baseDir]);
dirList = dirList([dirList.isdir]);
dirList = dirList(3:end);

for i = 1:length(dirList)
    subDirList = dir([baseDir dirList(i).name '\33*_RGB_RECT.t*']);
    subDirListComplete = dir([baseDir dirList(i).name]);
    subDirListComplete = subDirListComplete(3:end);
    if (length(subDirList)~=length(subDirListComplete))
        fprintf(dirList(i).name);
%         break;
    end
end

%%
fileName = 'D:\Data\Development\Projects\MSc GeoInformatics\Code\Results\MappingToolboxExp\MCD43A4.A2010033.h19v12.005.2010052115603.hdf';
info = hdfinfo(fileName);
hdftool(fileName)

band1 = hdfread('D:\Data\Development\Projects\MSc GeoInformatics\Code\Results\MappingToolboxExp\MCD43A4.A2010033.h19v12.005.2010052115603.hdf',...
    'MOD_Grid_BRDF', 'Fields', 'Nadir_Reflectance_Band1');
imagesc(band1)


%%
% Check out reverb hdf files
close all;clear all;
filename = 'G:\MSc GeoInformatics\Data\MODIS\MCD43A4.A2010025.h19v12.005.2010043064233.hdf';
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

info = hdfinfo(filename)
info.Attributes(2).Value

bands = [1 4 3 2];
figure;
clear im;
for i = 1:length(bands)
    data = hdfread(filename, 'MOD_Grid_BRDF', 'Fields', sprintf('Nadir_Reflectance_Band%d', bands(i)));
    nonodata = data(data~=(2^15-1));
    nonodata = double(nonodata)./(2^15-1);
    data(data==(2^15-1))=0;
    data = double(data)./(2^15-1);
    im(:,:,i) = imadjust(data, stretchlim(nonodata,[0 .99]), [0 1], 1);    
    subplot(2,2,i)
    imagesc(data);
end

rgbIm = im(:,:,1:3);
cirIm = im(:,:,[4 1 2]);

figure;
subplot(1,2,1)
imagesc(rgbIm)
subplot(1,2,2)
imagesc(cirIm)

R = maprasterref('RasterSize', [size(rgbIm, 1) size(rgbIm, 2)], ...
      'YLimWorld', [-4447802.078667 -3335851.559000 ], 'ColumnsStartFrom','north', ...
      'XLimWorld', [1111950.519667 2223901.039333])

tiffFileName = 'F:\MSc GeoInformatics\Data\MODIS\MCD43A4.A2010025.h19v12.005.2010043064233.tif';
geotiffwrite(tiffFileName, rgbIm, R, 'CoordRefSysCode', 'EPSG:27700'); %googled EPSG:4326 = sinusoid from

[rgbIm_ R_] =  geotiffread(tiffFileName);
gi = geotiffinfo(tiffFileName);
figure;
imagesc(rgbIm_)




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
close all;clear all;
modisFilename = 'F:\MSc GeoInformatics\Data\MODIS\MCD43A4.A2010025.h19v12.005.2010043064233.hdf';

% ngiFileName = 'G:\MSc GeoInformatics\Data\NGI\Rectified\RGB\3320AA\3320AA_01_2010_312_RGB_RECT.tif';
pciFileName = 'G:\MSc GeoInformatics\Data\NGI\My Rectified\3321D_2010_319\o3321D_2010_319_01_0004_RGBIR.tif';

% dos(sprintf('gdalinfo -proj4 "%s"', ngiFileName)) %'+proj=tmerc +lat_0=0 +lon_0=21 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs ' 
dos(sprintf('gdalinfo -proj4 "%s"', pciFileName)) %+proj=tmerc +lat_0=0 +lon_0=21 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs '

%Warp modis bands to tiff files in pci/ngi projection
bands = [1 4 3 2];
[p n] = fileparts(modisFilename);
for i = 1:length(bands)
    outFileName{i} = sprintf('%s/%s.Band%d.test.tif', p, n(1:end-4), bands(i));
    dos(sprintf('gdalwarp -r bilinear -t_srs "%s" HDF4_EOS:EOS_GRID:"%s":MOD_Grid_BRDF:Nadir_Reflectance_Band%d "%s"', ...        
        '+proj=tmerc +lat_0=0 +lon_0=21 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs', modisFilename, ...
        bands(i), outFileName{i}));
end

%Read tiff modis band images into one image

for i = 1:length(bands)
    [modisIm(:,:,(i)) r{i}] = geotiffread(outFileName{i});
end

pciInfo = geotiffinfo(pciFileName)
modisInfo = geotiffinfo(outFileName{1})

rgbnOutFileName = sprintf('%s/%s.RGBN.tif', p, n);

geotiffwrite(rgbnOutFileName, modisIm, r{end}, 'GeoKeyDirectoryTag', ...
    modisInfo.GeoTIFFTags.GeoKeyDirectoryTag);

%%
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
