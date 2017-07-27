close all; clear all;
cd 'D:\Data\Development\Projects\MSc GeoInformatics\Code\Results\BatchProcess'

imFileNames = {'D:\Data\Development\Projects\MSc GeoInformatics\Code\Results\BatchProcess\3321A_2010_316_01_0003_RGB2.tif', ...
    'D:\Data\Development\Projects\MSc GeoInformatics\Code\Results\BatchProcess\3321A_2010_316_01_0003_CIR2.tif'};
% imFileNames = { 'F:\MSc GeoInformatics\Data\NGI\Calibrated\3321C_2010_318\RGB\3321C_2010_318_01_0001_RGB.tif'};

% lutFileNames = {'F:\MSc GeoInformatics\Data\NGI\LUT\3321D\rgb_12.lut',...
%     'F:\MSc GeoInformatics\Data\NGI\LUT\3321C\rgb_12.lut'};

% icons = {'-x', '-o'};
for i = 1:length(imFileNames)
%     tmp = imread(imFileNames{i});
    t{i} = Tiff(imFileNames{i}, 'r+');
%     t{i}.setTag('BitsPerSample', 16);
    im{i} = t{i}.read();
end  

tout = Tiff('3321A_2010_316_01_0003_RGBIR_Deflate2.tif','w');

tagId = Tiff.TagID;
fn = fieldnames(tagId);
tagStruct = tagId;
for i = 1:length(fn)
    try
        val = t{1}.getTag(fn{i});
        tout.setTag(fn{i}, val);
        fn{i}, val
    catch
    end
end

imout = im{1};
imout = 16*cat(3, imout, im{2}(:,:,1));
tout.setTag('BitsPerSample', 16);
tout.setTag('SamplesPerPixel', 4);
tout.setTag('MaxSampleValue', 65535);
% tout.setTag('MaxSampleValue', 4095);
% tout.setTag('NoDataValue', 4095);
tout.setTag('Compression', Tiff.Compression.Deflate); %DEFLATE IS THE PNG ALGORITHM
% tout.rewriteDirectory();
tout.write(imout);
% tout.setTag('BitsPerSample', 12);
% tout.rewriteDirectory();
tout.close();
t{1}.close();
t{2}.close();

%Check if CIR+B from RGB works ok

rg/rg2
gg/gg2

im1 = im{1}(1:2000,1:2000,:);
im2 = im{2}(1:2000,1:2000,:);

rg = median(mean(im1(:,:,1)))
gg = median(mean(im1(:,:,2)))

rg2 = median(mean(im2(:,:,2)))
gg2 = median(mean(im2(:,:,3)))

im4band = cat(3, (rg/rg2)*im2(:,:,2), (gg/gg2)*im2(:,:,3),  im1(:,:,3), im2(:,:,1));

figure;
h(1) = subplot(1,3,1);
imshow(16*im4band(:,:,[1:3]));
h(2) = subplot(1,3,2);
imshow(16*im4band(:,:,[4 1 2]));
h(3) = subplot(1,3,3);
imshow(16*im1);
linkaxes(h, 'xy')

%%
% Start the real thing
close all; clear all;
inBaseDirName = 'F:\MSc GeoInformatics\Data\NGI\Calibrated\3321D_2010_319';
outDirName = 'G:\MSc GeoInformatics\Data\NGI\My Calibrated\3321D_2010_319';

if (~exist(outDirName, 'dir'))
    mkdir(outDirName);
end
BatchMake4BandTiff(inBaseDirName, outDirName);

% im4band = cat(3, im4band, );
% tout.write(imout);
%NOTES
%--------------------------------------------------------------------------
%- 12bit jpegs cannot be loaded by the default libjpeg which is linked to
%libtiff.  A special 12 bit libjpeg needs to be compiled (this will then
%only support 12bit jpegs and not 8bit).  libtiff can then be linked with
%the 12 bit version.  GDAL apparently is now linked with both 8 and 12bit 
%versions of libjpeg. 
%- We have a few options for dealing with 12bit tiffs in matlab: 
%(1) convert them outside of matlab using gdal / arcmap 
%(2) get a 12bit version of libtiff and replace the one Matlab is using
%with it (its a mex file so we'd have to recompile ourselves).  Bad.
%(3) write gdal code to do everything we want to do in matlab i.e.
%make 4 band, 16bit tiff
%- Different RG LUT's are used for CIR and RGB
%-NB the real res is better than 0.5m so rectifying to 0.5m is actually
%downsampling a bit
%-NB packing the 12bits into 16bit channels but with BitsPerSample set to
%12 gives the correct max val = 4096 in arcmap and thus renders "out the
%box"

% gdal_translate -ot UInt16 -of GTiff -co "TILED=YES" -co "COMPRESS=DEFLATE" 3321A_2010_316_01_0003_RGB.tif 3321A_2010_316_01_0003_RGB_16b.tif